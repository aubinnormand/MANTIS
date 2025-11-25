import os
import math
import numpy as np
import geopandas as gpd
from shapely.geometry import box
import matplotlib.pyplot as plt
from pyproj import Geod

# ------------------------- #
#         UTILITAIRES       #
# ------------------------- #
# initialiser l'ellipsoïde WGS84
geod = Geod(ellps="WGS84")
def geod_area(poly):
    # calcule l'aire en m² (la fonction renvoie aussi un périmètre, que l'on ignore)
    area, _ = geod.geometry_area_perimeter(poly)
    return abs(area)  # aire positive
        
def simplify_world_geometry(world_gdf, tolerance=0.1):
    """
    Simplifie la géométrie du GeoDataFrame mondial pour accélérer les calculs et l'affichage.
    
    - world_gdf : GeoDataFrame du monde
    - tolerance : tolérance en degrés (~0.1° correspond à environ 11 km)
    """
    simplified = world_gdf.copy()
    simplified["geometry"] = simplified["geometry"].simplify(tolerance, preserve_topology=True)
    return simplified

def compute_degrees_per_km(latitude):
    """Retourne le nombre de degrés en latitude et longitude pour 1 km."""
    lat_deg = 1 / 111.32
    lon_deg = 1 / (111.32 * math.cos(math.radians(latitude)))
    return lat_deg, lon_deg

def ensure_unique_spatial_keys(grid, cle_geo):
    """
    Rend les clés uniques en ajoutant un suffixe si nécessaire.
    Corrige le FutureWarning de Pandas.
    """
    grid["tmp"] = grid.groupby(cle_geo).cumcount()
    
    # Convertir la colonne en str avant concaténation
    grid[cle_geo] = grid[cle_geo].astype(str)
    mask = grid["tmp"] > 0
    grid.loc[mask, cle_geo] = grid.loc[mask, cle_geo] + "_" + grid.loc[mask, "tmp"].astype(str)
    
    grid.drop(columns="tmp", inplace=True)

def report_duplicate_keys(grid, cle_geo):
    """Affiche le nombre de doublons dans la clé géographique."""
    dup = grid[cle_geo].duplicated().sum()
    print(f"🔍 {dup} doublon(s)" if dup else "✔ Aucun doublon")


# ------------------------- #
#     CREATION DE GRILLE    #
# ------------------------- #

import numpy as np
import geopandas as gpd
from shapely.geometry import box

def build_country_grid_WGS84(gdf_country, grid_size_km, cle_geo):
    """
    Construit une grille régulière en WGS84 sur un pays donné.

    Paramètres
    ----------
    gdf_country : GeoDataFrame
        GeoDataFrame du pays avec géométrie.
    grid_size_km : float
        Taille de la maille en kilomètres.
    cle_geo : str
        Nom de la colonne pour la clé unique de maille.

    Retour
    ------
    GeoDataFrame
        Grille découpée par le pays avec cle_geo, rank_x, rank_y.
    """
    # Union des géométries du pays
    geom = gdf_country.geometry.unary_union
    minx, miny, maxx, maxy = geom.bounds

    # Calcul degrés/km pour la latitude moyenne
    midpoint_lat = (miny + maxy) / 2
    lat_deg, lon_deg = compute_degrees_per_km(midpoint_lat)
    dx, dy = grid_size_km * lon_deg, grid_size_km * lat_deg

    # Alignement sur 0° lat/lon
    minx_aligned = np.floor(minx / dx) * dx
    miny_aligned = np.floor(miny / dy) * dy

    # Création des coordonnées
    xs = np.arange(minx_aligned, maxx + dx, dx)
    ys = np.arange(miny_aligned, maxy + dy, dy)

    # Création des polygones et clés
    rects = []
    cle_list = []
    rank_x_list = []
    rank_y_list = []

    for i, x in enumerate(xs[:-1]):
        for j, y in enumerate(ys[:-1]):
            rects.append(box(x, y, x + dx, y + dy))
            cle_list.append(f"{i}_{j}")
            rank_x_list.append(i)
            rank_y_list.append(j)

    # Construction du GeoDataFrame
    grid = gpd.GeoDataFrame(
        {
            cle_geo: cle_list,
            'rank_x': rank_x_list,
            'rank_y': rank_y_list,
            'geometry': rects
        },
        crs=gdf_country.crs
    )

    # Intersection avec le pays
    grid = gpd.overlay(grid, gdf_country, how="intersection")

    # Supprimer les doublons au cas où l'overlay crée plusieurs morceaux
    grid = grid.drop_duplicates(subset=[cle_geo]).reset_index(drop=True)

    return grid


# ------------------------- #
#   FONCTION PRINCIPALE     #
# ------------------------- #
def generate_country_grid(
        world_gdf,
        country_name,
        name_attribute,
        grid_size_km,
        cle_geo,
        code_attribute=None):
    """
    Fonction unique à appeler :
    - sélectionne un pays
    - génère sa grille WGS84
    - rend les clés uniques
    - contrôle les doublons
    - préfixe la clé géographique par le code de zone
    - sauvegarde en GeoJSON
    """

    # 1) Sélection du pays
    country = world_gdf[world_gdf[name_attribute] == country_name]
    if country.empty:
        raise ValueError(f"❌ Pays '{country_name}' non trouvé dans '{name_attribute}'")

    # 2) Détermination du code de zone
    if code_attribute and code_attribute in country.columns:
        zone_code = str(country.iloc[0][code_attribute])
    else:
        zone_code = input(f"Entrer le code de zone pour {country_name} : ").strip()

    # 3) Construction
    grid = build_country_grid_WGS84(country, grid_size_km, cle_geo)
    
     # 4) Préfixer la clé par le code de zone 
    grid[cle_geo] = zone_code + "_" + grid[cle_geo].astype(str)
    
    # 5) Calculer l'aire des cellules
    grid=compute_area_grid(grid)

    # 5) Correction clés pour doublons éventuels
    #ensure_unique_spatial_keys(grid, cle_geo)

    # 6) Vérification
    report_duplicate_keys(grid, cle_geo)

    # 7) Conserver seulement cle_geo et geometry
    grid = grid[[cle_geo, 'geometry','area_km2','rank_x','rank_y']]

    return grid

def save_country_grid(grid_gdf, path_output, country_name, cle_geo):
    filename = f"{country_name}_{cle_geo}.geojson"
    
    """Crée un dossier si nécessaire et exporte la grille en GeoJSON."""
    os.makedirs(path_output, exist_ok=True)
    
    grid_gdf.to_file(os.path.join(path_output, filename), driver="GeoJSON")

    print(f"📍 Grille enregistrée : {os.path.join(path_output, filename)}")

def plot_country_grid(world_gdf, grid_gdf, name_attribute, country_name):
    import matplotlib.pyplot as plt

    country = world_gdf[world_gdf[name_attribute] == country_name]
    if country.empty:
        raise ValueError(f"❌ Pays '{country_name}' non trouvé")

    fig, ax = plt.subplots(figsize=(10, 8))
    world_gdf.plot(ax=ax, color="lightgrey", edgecolor="white", linewidth=0.3)
    country.plot(ax=ax, color="#69b3e7", edgecolor="black", linewidth=1)
    grid_gdf.boundary.plot(ax=ax, color="red", linewidth=0.5)

    minx, miny, maxx, maxy = grid_gdf.total_bounds
    ax.set_xlim(minx, maxx)
    ax.set_ylim(miny, maxy)
    ax.set_title(f"Grille sur {country_name}")
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")

    plt.tight_layout()
    plt.show() 


def compute_area_grid(grid):
    # initialiser l'ellipsoïde WGS84 (le même qu’EPSG:4326)
    geod = Geod(ellps="WGS84")
    grid = grid.to_crs(epsg=4326)
    
    # appliquer directement sur ton dataframe "carte_maille"
    grid["area_km2"] = grid.geometry.apply(geod_area)/ 1e6
    return grid
