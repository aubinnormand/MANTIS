import os
import numpy as np
import pandas as pd
import geopandas as gpd

# -------------------------
# Paramètres globaux
# -------------------------
CLE_GEO_DEFAULT = "codeMaille10Km"
CLE_ID_DEFAULT = "speciesKey"
PATH_DATA = r"C:\path\to\data"  # à adapter
SOURCE = "GBIF"
VAR_OBS_DEFAULT = "nombreObs"
SEUIL_FUSION_DEFAULT = 1000
METHODE_FUSION_DEFAULT = "barycentre"
MAX_DISTANCE_DEFAULT = 100

# -------------------------
# Fonctions de fusion de cellules
# -------------------------

def fusion_cells_by_obs(
        gdf,
        cle_geo=CLE_GEO_DEFAULT,
        var_obs=VAR_OBS_DEFAULT,
        seuil=SEUIL_FUSION_DEFAULT,
        methode=METHODE_FUSION_DEFAULT,
        max_distance=MAX_DISTANCE_DEFAULT,
        display=True):
    """
    Fusionne les cellules tant que var_obs < seuil et garde la trace des cellules fusionnées.
    """
    gdf = gdf.copy()
    n_fusion = 0
    gdf['fusion'] = gdf[cle_geo].apply(lambda x: [x])

    if methode == "barycentre":
        gdf['centroid'] = gdf['geometry'].centroid

    while len(gdf) > 1 and gdf[var_obs].min() < seuil:
        n_fusion += 1
        if display:
            print(f"Fusion n°{n_fusion}, n_min= {gdf[var_obs].min()}")
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

        # Fusion
        gdf.at[neighbor_idx, var_obs] += gdf.at[idx_min, var_obs]
        gdf.at[neighbor_idx, 'geometry'] = gdf.at[neighbor_idx, 'geometry'].union(gdf.at[idx_min, 'geometry'])
        gdf.at[neighbor_idx, 'fusion'] += gdf.at[idx_min, 'fusion']
        gdf = gdf.drop(idx_min)

        if methode == "barycentre":
            gdf.loc[neighbor_idx, 'centroid'] = gdf.loc[neighbor_idx, 'geometry'].centroid

    if 'centroid' in gdf.columns:
        gdf = gdf.drop(columns=['centroid'])
    return gdf.reset_index(drop=True)


def apply_fusion_to_biodiv(
        df_biodiv,
        merged_fusion,
        cle_geo=CLE_GEO_DEFAULT,
        cle_ID=CLE_ID_DEFAULT,
        var_obs=VAR_OBS_DEFAULT):
    """
    Applique les fusions de cellules à df_biodiv et regroupe par cellule, espèce et période.
    """
    df_dico = generer_dictionnaire_taxonomie(df_biodiv, cle_ID)

    mapping = {}
    for _, row in merged_fusion.iterrows():
        final_cell = row[cle_geo]
        for cell in row['fusion']:
            mapping[cell] = final_cell

    mapping_df = pd.DataFrame(list(mapping.items()), columns=['orig', 'final'])
    df_biodiv[cle_geo] = df_biodiv[cle_geo].astype(str)
    mapping_df['orig'] = mapping_df['orig'].astype(str)
    
    df_merged = df_biodiv.merge(mapping_df, left_on=cle_geo, right_on='orig', how='left')
    df_merged[cle_geo] = df_merged['final']
    df_merged = df_merged.drop(columns=['orig', 'final'])
    
    df_grouped = df_merged.groupby([cle_geo, cle_ID, 'periode'], as_index=False)[var_obs].sum()
    df_grouped = pd.merge(df_grouped, df_dico, on=cle_ID)
    return df_grouped


def fusionner_depuis_fichiers(
        country_name,
        grid_size_km,
        domaine,
        bornes_temporelles,
        path_data=PATH_DATA,
        cle_geo=CLE_GEO_DEFAULT,
        cle_ID=CLE_ID_DEFAULT,
        source=SOURCE,
        var_obs=VAR_OBS_DEFAULT,
        seuil_fusion=SEUIL_FUSION_DEFAULT,
        methode_fusion=METHODE_FUSION_DEFAULT,
        max_distance=MAX_DISTANCE_DEFAULT):
    """
    Effectue la fusion des mailles à partir des fichiers déjà traités pour un domaine donné.
    """
    # Charger la grille
    grid_path = generate_grid_path(country_name, grid_size_km, domaine, cle_geo, source=source, fusion=False)
    grid = gpd.read_file(grid_path)

    # Charger df_final
    chaine_bornes = "_".join(map(str, bornes_temporelles))
    initial_file_path = os.path.join(
        path_data, source, "processed", f"{source}_{country_name}",
        f"data_{source}_{country_name}_{domaine}_{cle_geo}_{cle_ID}_periodes{chaine_bornes}.csv"
    )
    df_final = pd.read_csv(initial_file_path)

    # Somme des observations par cellule
    df_obs = df_final.groupby(cle_geo, as_index=False)[var_obs].sum()
    merged = grid.merge(df_obs, on=cle_geo, how="left")
    merged[var_obs] = merged[var_obs].fillna(0).astype(float)

    # Appliquer la fusion
    merged_fusion = fusion_cells_by_obs(
        merged, cle_geo=cle_geo, var_obs=var_obs, seuil=seuil_fusion,
        methode=methode_fusion, max_distance=max_distance
    )

    # Appliquer les résultats à df_biodiv
    df_grouped = apply_fusion_to_biodiv(df_final, merged_fusion, cle_geo=cle_geo, cle_ID=cle_ID, var_obs=var_obs)
    return df_grouped
