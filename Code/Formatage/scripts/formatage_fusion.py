import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
import time
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from formatage_ObsToGrid import generate_taxo_dict

def fusion_cells_by_obs(df_final, grid, cle_geo, var_obs, seuil=1000, methode="barycentre", max_distance=100, display=True):
    """
    Fusionne les cellules tant que nombreObs < seuil et garde la trace des cellules fusionnées.
    """
    # 1️⃣ somme des observations par cellule
    df_obs = df_final.groupby(cle_geo, as_index=False)[var_obs].sum()
    
    # 2️⃣ merge avec la grille
    gdf = grid.merge(df_obs, on=cle_geo, how="left")
    gdf[var_obs] = gdf[var_obs].fillna(0).astype(float)
    
    # colonne pour garder la trace des cellules fusionnées
    gdf['fusion'] = gdf[cle_geo].apply(lambda x: [x])
    
    # barycentre si nécessaire
    if methode == "barycentre":
        gdf['centroid'] = gdf['geometry'].centroid

    n_init = len(gdf)
    n_last_display = n_init

    print("🔄 Fusion des cellules en cours...") 
    
    while len(gdf) > 1 and gdf[var_obs].min() < seuil:
        idx_min = gdf[var_obs].idxmin()
        cell = gdf.loc[idx_min]
    
        # --- logique de fusion (barycentre ou autre) ---
        # (reste identique à ce que vous avez)
    
        # après fusion et suppression de la cellule
        n_current = len(gdf)
        pct_remaining = n_current / n_init
    
        # afficher seulement si la réduction >=1% par rapport au dernier affichage
        if display and n_current <= n_last_display * 0.90:
            print(f"{n_current} cellules restantes ({pct_remaining:.1%}), min {var_obs} = {gdf[var_obs].min()}")
            n_last_display = n_current

        idx_min = gdf[var_obs].idxmin()
        cell = gdf.loc[idx_min]

        if methode == "barycentre":
            centroid_cell = cell['centroid']
            others = gdf.drop(idx_min)
            if others.empty:
                break
            others['distance'] = others['centroid'].apply(lambda c: centroid_cell.distance(c))
            neighbor_idx = others['distance'].idxmin()
        else:
            x_min, y_min = cell['rank_x'], cell['rank_y']
            neighbors = pd.DataFrame()
            for dist in range(1, max_distance+1):
                neighbors = gdf[
                    ((abs(gdf['rank_x'] - x_min) == dist) & (gdf['rank_y'] == y_min)) |
                    ((abs(gdf['rank_y'] - y_min) == dist) & (gdf['rank_x'] == x_min))
                ].drop(idx_min, errors='ignore')
                neighbors = neighbors[neighbors[var_obs].notna()]
                if not neighbors.empty:
                    break
            if neighbors.empty:
                break

            if methode == "max":
                val = neighbors[var_obs].max()
                candidates = neighbors[neighbors[var_obs] == val]
            elif methode == "min":
                val = neighbors[var_obs].min()
                candidates = neighbors[neighbors[var_obs] == val]
            elif methode == "random":
                candidates = neighbors
            else:
                raise ValueError("methode doit être 'max', 'min', 'random' ou 'barycentre'")
            
            neighbor_idx = np.random.choice(candidates.index)

        # Fusionner
        gdf.at[neighbor_idx, var_obs] += gdf.at[idx_min, var_obs]
        gdf.at[neighbor_idx, 'geometry'] = gdf.at[neighbor_idx, 'geometry'].union(gdf.at[idx_min, 'geometry'])
        gdf.at[neighbor_idx, 'fusion'] += gdf.at[idx_min, 'fusion']
        gdf = gdf.drop(idx_min)

        if methode == "barycentre":
            gdf.loc[neighbor_idx, 'centroid'] = gdf.loc[neighbor_idx, 'geometry'].centroid

    if 'centroid' in gdf.columns:
        gdf = gdf.drop(columns=['centroid'])
        
    print(f"✅ Fusion terminée ! Nombre final de cellules : {len(gdf)}")  # Fin de la fusion

    return gdf.reset_index(drop=True)

def save_merged_grid(merged_fusion, path_data, zone, cle_geo, source, seuil_fusion):
    
    # Colonnes à garder
    colonnes_a_garder = [cle_geo, "geometry",'area_km2']
    
    grid_fusion = merged_fusion[colonnes_a_garder].copy()

    # Génération du chemin du fichier
    filename = f"{zone}_{cle_geo}_fused_{source}_{seuil_fusion}.geojson"
    path_output = path_data / "SIG_zone" / zone /filename

    # Sauvegarde de la grille fusionnée
    """Crée un dossier si nécessaire et exporte la grille en GeoJSON."""
    path_output.parent.mkdir(parents=True, exist_ok=True)
    
    grid_fusion.to_file(path_output, driver="GeoJSON")

    print(f"📍 Grille enregistrée : {path_output}")
    
    return path_output



def apply_fusion_to_biodiv(df_biodiv, merged_fusion, cle_geo='codeMaille10Km', cle_ID='speciesKey', var_obs='nombreObs'):
    """
    Applique les fusions à df_biodiv et sauvegarde le résultat si path_save est fourni.
    """
    df_dico = generate_taxo_dict(df_biodiv, cle_ID)

    # Mapping original -> final
    mapping = {cell: row[cle_geo] for _, row in merged_fusion.iterrows() for cell in row['fusion']}
    mapping_df = pd.DataFrame(list(mapping.items()), columns=['orig', 'final'])

    df_biodiv[cle_geo] = df_biodiv[cle_geo].astype(str)
    mapping_df['orig'] = mapping_df['orig'].astype(str)
    df_merged = df_biodiv.merge(mapping_df, left_on=cle_geo, right_on='orig', how='left')
    df_merged[cle_geo] = df_merged['final']
    df_merged = df_merged.drop(columns=['orig', 'final'])

    # Regrouper
    df_grouped = df_merged.groupby([cle_geo, cle_ID, 'year', 'month'], as_index=False)[var_obs].sum()
    df_grouped = pd.merge(df_grouped, df_dico, on=cle_ID)

    return df_grouped

def plot_fusionned_grid(gdf, var_obs='nombreObs', title=None, cmap='viridis', figsize=(10,10), log_scale=False):
    """
    Affiche une grille fusionnée avec un dégradé de couleur représentant le nombre d'observations.
    
    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame fusionné contenant la colonne `var_obs`.
    var_obs : str
        Nom de la colonne contenant le nombre d'observations.
    title : str
        Titre du graphique.
    cmap : str
        Palette de couleur pour le plotting.
    figsize : tuple
        Taille de la figure.
    log_scale : bool
        Si True, utilise une échelle logarithmique pour les couleurs.
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Choisir normalisation
    if log_scale:
        norm = colors.LogNorm(vmin=max(gdf[var_obs].min(), 1), vmax=gdf[var_obs].max())
    else:
        norm = None

    gdf.plot(column=var_obs, ax=ax, cmap=cmap, legend=True, norm=norm,
             edgecolor='black', linewidth=0.5)

    ax.set_axis_off()
    if title:
        ax.set_title(title, fontsize=16)
    plt.show()


def save_fusionned_biodiv(df_final, source, zone, cle_geo,path_data,seuil):
    output_path = path_data / source / 'processed' / zone/f"{source}_{zone}_{cle_geo}_fused{str(seuil)}.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_final.to_csv(output_path, index=False)
    print(f"🎉 Données nettoyées sauvegardées dans : {output_path}")
    return output_path
