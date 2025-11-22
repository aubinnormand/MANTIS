import os
import math
import numpy as np
import geopandas as gpd
from shapely.geometry import box
import matplotlib.pyplot as plt

# -------------------------
# Paramètres globaux
# -------------------------
CLE_GEO_DEFAULT = "codeMaille100Km"
WORLD_TERRESTRE_FILE = "geoBoundariesCGAZ_ADM0.shp"
WORLD_MARITIME_FILE = "eez_v11.gpkg"
CRITERE_MARITIME_DEFAULT = "ISO_TER1"
CRITERE_TERRESTRE_DEFAULT = "shapeGroup"
PATH_DATA = r"C:\path\to\data"  # à adapter
SOURCE = "GBIF"


# -------------------------
# Fonctions d'importation
# -------------------------
def charger_donnees_geo(base_path=PATH_DATA, source=SOURCE,
                        terrestrial_file=WORLD_TERRESTRE_FILE,
                        maritime_file=WORLD_MARITIME_FILE):
    """
    Charge les données géographiques et retourne les GeoDataFrames pertinents.
    """
    geo_data = load_geospatial_data(base_path, source, terrestrial_file, maritime_file)
    return geo_data["world_terrestre"], geo_data["world_maritime"]

    
def rename_sig_files(country_name, grid_size_km, path_data=PATH_DATA, source=SOURCE):
    """
    Duplique et renomme les fichiers SIG en ajoutant le suffixe "_rec".
    """
    sig_folder = os.path.join(path_data, source, "processed", f"{source}_{country_name}")

    if not os.path.exists(sig_folder):
        print(f"⚠️ Le dossier SIG n'existe pas : {sig_folder}")
        return
        
def load_geospatial_data(path_sig, source="GBIF",
                         terrestrial_file=WORLD_TERRESTRE_FILE,
                         maritime_file=WORLD_MARITIME_FILE):
    """
    Charge les fichiers SIG globaux (terrestre et maritime).

    Args:
        base_path (str): Chemin du répertoire du projet.
        source (str): Source des données, par défaut "GBIF".
        terrestrial_file (str): Nom du fichier terrestre.
        maritime_file (str): Nom du fichier maritime.

    Returns:
        dict: Un dictionnaire contenant les GeoDataFrames chargés.
    """
    geodata = {
        "world_terrestre": gpd.read_file(os.path.join(path_sig, terrestrial_file)),
        "world_maritime": gpd.read_file(os.path.join(path_sig, maritime_file)),
    }

    return geodata
# -------------------------
# Fonctions utilitaires
# -------------------------
def degrees_per_km(latitude):
    """Calcule la conversion degrés/km pour une latitude donnée."""
    lat_deg_per_km = 1 / 111.32
    lon_deg_per_km = 1 / (111.32 * math.cos(math.radians(latitude)))
    return lat_deg_per_km, lon_deg_per_km


def make_geo_keys_unique(grid, cle_geo=CLE_GEO_DEFAULT):
    """Rend les clés de mailles uniques en ajoutant un suffixe si nécessaire."""
    if not grid.empty:
        grid["cle_geo_unique"] = grid.groupby(cle_geo).cumcount().astype(str)
        grid.loc[grid["cle_geo_unique"] != "0", cle_geo] += "_" + grid["cle_geo_unique"]
        grid.drop(columns=["cle_geo_unique"], inplace=True)


def check_duplicates(grid, grid_name, cle_geo=CLE_GEO_DEFAULT):
    """Vérifie les doublons pour une clé de grille donnée et affiche un message."""
    if not grid.empty:
        duplicate_count = grid[cle_geo].duplicated().sum()
        if duplicate_count > 0:
            print(f"⚠️ Attention : {duplicate_count} doublon(s) trouvé(s) dans {grid_name} !")
        else:
            print(f"✅ Aucun doublon trouvé dans {grid_name}.")

# -------------------------
# Fonctions de création de grille
# -------------------------
def creer_grille_pays(
    path_data,
    country_name,
    grid_size_km,
    world_terrestre,
    world_maritime,
    cle_geo=CLE_GEO_DEFAULT,
    name_attribute="shapeName",
    code_attribute="shapeGroup",
    source=SOURCE
):
    """
    Crée la grille pour un pays et charge les fichiers GeoJSON correspondants.
    """
    country_code = world_terrestre[world_terrestre[name_attribute] == country_name][code_attribute].iloc[0]

    # Génération des grilles
    generer_grille_pays(
        path_data, country_name, grid_size_km,
        world_terrestre, world_maritime,
        critere_terrestre=code_attribute,
        critere_maritime=CRITERE_MARITIME_DEFAULT,
        name_attribute="shapeName",
        code_attribute="shapeGroup",
        source=source
    )

    # Chargement des fichiers générés
    path_grid_terrestre = generate_grid_path( country_name, grid_size_km, 'terrestre', cle_geo=cle_geo,path_data=path_data, source=source)
    path_grid_maritime = generate_grid_path( country_name, grid_size_km, 'maritime', cle_geo=cle_geo, path_data=path_data,source=source)
    path_grid_combined = generate_grid_path(country_name, grid_size_km, 'combined', cle_geo=cle_geo, path_data=path_data,source=source)

    country_grid_terrestre = load_grid_file(path_grid_terrestre, "Grid Terrestre")
    country_grid_maritime = load_grid_file(path_grid_maritime, "Grid Maritime")
    country_grid_combined = load_grid_file(path_grid_combined, "Grid Combined")

    return country_grid_terrestre, country_grid_maritime, country_grid_combined
    
def generer_grille_pays(
    path_data,
    country_name,
    grid_size_km,
    world_terrestre,
    world_maritime,
    critere_terrestre=CRITERE_TERRESTRE_DEFAULT,
    critere_maritime=CRITERE_MARITIME_DEFAULT,
    name_attribute="shapeName",
    code_attribute="shapeGroup",
    source=SOURCE
):
    """
    Génère les grilles terrestre, maritime et combinée pour un pays donné,
    puis les sauvegarde en GeoJSON.
    """
    cle_geo = f"codeMaille{grid_size_km}Km"

    # 1. Identifier le code du pays
    country_code = world_terrestre[world_terrestre[name_attribute] == country_name][code_attribute].iloc[0]

    # 2. Sélection des géométries terrestre et maritime
    country_terrestre = world_terrestre[world_terrestre[critere_terrestre] == country_code]
    country_maritime = world_maritime[world_maritime[critere_maritime] == country_code]

    # Fusion des géométries
    geom_terrestre = country_terrestre.geometry.unary_union
    geom_maritime = country_maritime.geometry.unary_union if not country_maritime.empty else None
    geom_fusionnee = geom_terrestre.union(geom_maritime) if geom_maritime is not None else geom_terrestre

    # 3. Création des grilles
    country_grid_terrestre = create_and_rename_grid(
        country_terrestre, country_code, critere_terrestre, grid_size_km, cle_geo=cle_geo, crop=True, display=False
    )
    
    if cle_geo not in country_grid_terrestre.columns:
        print("⚠ Attention, cle_geo non trouvée :", cle_geo)
        print("Colonnes disponibles generer_grille_pays :", country_grid_terrestre.columns)

    country_grid_maritime = gpd.GeoDataFrame()
    if not country_maritime.empty:
        country_grid_maritime = create_and_rename_grid(
            country_maritime, country_code, critere_maritime, grid_size_km,cle_geo= cle_geo, crop=True, display=False
        )

    # Grille combinée
    gdf_fusionne = gpd.GeoDataFrame(geometry=[geom_fusionnee], crs=world_terrestre.crs)
    gdf_fusionne["Code"] = country_code
    country_grid_combined = create_and_rename_grid(
        gdf_fusionne, country_code, "Code", grid_size_km, cle_geo=cle_geo, crop=False, display=False
    )


    # Rendre les noms de mailles uniques
    make_geo_keys_unique(country_grid_terrestre, cle_geo=cle_geo)
    make_geo_keys_unique(country_grid_maritime,cle_geo= cle_geo)
    make_geo_keys_unique(country_grid_combined,cle_geo= cle_geo)

    # Vérification des doublons
    check_duplicates(country_grid_terrestre, "country_grid_terrestre",cle_geo= cle_geo)
    check_duplicates(country_grid_maritime, "country_grid_maritime", cle_geo=cle_geo)
    check_duplicates(country_grid_combined, "country_grid_combined", cle_geo=cle_geo)

    # 4. Sauvegarde des grilles
    country_grid_terrestre.to_file(
        generate_grid_path( country_name, grid_size_km, 'terrestre', cle_geo=cle_geo, path_data=path_data,source=source),
        driver="GeoJSON"
    )

    if not country_grid_maritime.empty:
        country_grid_maritime.to_file(
            generate_grid_path(country_name, grid_size_km, 'maritime', cle_geo=cle_geo,path_data=path_data,source=source),
            driver="GeoJSON"
        )

    country_grid_combined.to_file(
        generate_grid_path(country_name, grid_size_km, 'combined', cle_geo=cle_geo, path_data=path_data,source=source),
        driver="GeoJSON"
    )
def create_and_rename_grid(source_gdf, country_code, column_name, grid_size_km, cle_geo=CLE_GEO_DEFAULT, crop=True, display=False):

    # Calcul du centre de latitude pour ajuster la grille si nécessaire
    if source_gdf.geometry.is_empty.all():
        midpoint_lat = None
    else:
        minx, miny, maxx, maxy = source_gdf.geometry.unary_union.bounds
        midpoint_lat = (miny + maxy) / 2

    # Création de la grille
    grid = create_country_grid_WGS84(
        source_gdf,
        country_code,
        column_name,
        grid_size_km,
        midpoint_lat=midpoint_lat,
        crop=crop,
        display=display
    )

    print("Colonnes disponibles create_and_rename_grid:", grid.columns)

    # Renommage de la colonne
    grid = grid.rename(columns={'cell_name': cle_geo})

    return grid

    
def create_country_grid_WGS84(
        gdf,
        country_code,
        col_code=CRITERE_TERRESTRE_DEFAULT,
        grid_size_km=100,
        midpoint_lat=None,
        crop=False,
        display=False,
        cle_geo=CLE_GEO_DEFAULT):
    """
    Génère une grille pour un pays spécifique en WGS84.
    """
    country = gdf[gdf[col_code] == country_code]
    if country.empty:
        raise ValueError(f"Pays '{country_code}' introuvable dans le GeoDataFrame.")

    def create_grid(country, grid_size, midpoint_lat=None, crop=False):
        """Création d'une grille avec une taille définie."""
        country_geometry = country.geometry.unary_union
        minx, miny, maxx, maxy = country_geometry.bounds
        
        if midpoint_lat is None:
            midpoint_lat = (miny + maxy) / 2

        lat_deg_per_km, lon_deg_per_km = degrees_per_km(midpoint_lat)
        dy, dx = grid_size * lat_deg_per_km, grid_size * lon_deg_per_km

        x_coords = np.arange(minx, maxx + dx, dx)
        y_coords = np.arange(miny, maxy + dy, dy)
        polygons = []
        for x0 in x_coords[:-1]:
            for y0 in y_coords[:-1]:
                polygons.append(box(x0, y0, x0 + dx, y0 + dy))

        grid_gdf = gpd.GeoDataFrame({cle_geo: range(len(polygons))}, geometry=polygons, crs=gdf.crs)

        if crop:
            grid_gdf = gpd.overlay(grid_gdf, country, how='intersection')

        if display:
            grid_gdf.boundary.plot()
            country.boundary.plot(edgecolor='red', ax=plt.gca())
            plt.show()
            
        print("Colonnes disponibles create_grid:", grid_gdf.columns)
        return grid_gdf
    grid=create_grid(country, grid_size_km, midpoint_lat, crop)
    print("Colonnes disponibles create_country_grid_WGS84:", grid.columns)

    return create_grid(country, grid_size_km, midpoint_lat, crop)


def calcul_largeur_maille_pays(country_name, world_terrestre, path_fichier, n_donnees_par_maille):
    """
    Calcule la largeur de maille recommandée pour un pays donné à partir d'un fichier CSV.
    """
    country_shape = world_terrestre[world_terrestre['shapeName'] == country_name]
    if country_shape.empty:
        raise ValueError(f"❌ Pays {country_name} introuvable dans world_terrestre")

    country_shape_metric = country_shape.to_crs(epsg=3857)
    surface_km2 = country_shape_metric.geometry.area.sum() / 10**6

    with open(path_fichier, "r", encoding="utf-8") as f:
        nb_lignes_sum = sum(1 for _ in f)

    densite_donnees = nb_lignes_sum / surface_km2
    largeur_maille_recommande = round(math.sqrt(n_donnees_par_maille / densite_donnees))
    print(f"La largeur de maille recommandée pour {country_name} est de {largeur_maille_recommande:.2f} km")
    return largeur_maille_recommande


def generate_grid_path(country_name, grid_size_km, grid_type, cle_geo=CLE_GEO_DEFAULT,path_data=PATH_DATA,
                       source=SOURCE, fusion=False, seuil_fusion=1000):
    """
    Génère le chemin complet pour un fichier de grille SIG.
    """
    print(country_name)
    sig_folder = os.path.join(path_data, source, "processed", f"{source}_{country_name}")
    print(sig_folder)
    os.makedirs(sig_folder, exist_ok=True)

    base_filename = f"grid_{country_name}_{grid_type}_{cle_geo}"
    filename = f"{base_filename}.geojson"
    if fusion:
        filename = f"{base_filename}Adapt{seuil_fusion}.geojson"

    return os.path.join(sig_folder, filename)


def load_grid_file(file_path, name="Grille"):
    """
    Vérifie si un fichier de grille existe et le charge.
    """
    if os.path.exists(file_path):
        print(f"✅ {name} trouvée et chargée avec succès !")
        return gpd.read_file(file_path)
    else:
        print(f"⚠️ {name} non trouvée. Elle ne sera pas traitée.")
        return None


def duplicate_grid_files(sig_folder, grid_size_km):
    """
    Duplique les fichiers SIG existants en ajoutant le suffixe "_rec".
    """
    for file_name in os.listdir(sig_folder):
        suffix = f"codeMaille{grid_size_km}Km.geojson"
        if file_name.endswith(suffix):
            old_path = os.path.join(sig_folder, file_name)
            new_file_name = file_name.replace(suffix, f"codeMaille{grid_size_km}Km_rec.geojson")
            new_path = os.path.join(sig_folder, new_file_name)

            with open(old_path, 'rb') as f_src:
                content = f_src.read()
            with open(new_path, 'wb') as f_dst:
                f_dst.write(content)

            print(f"✅ Fichier dupliqué : {file_name} → {new_file_name}")
