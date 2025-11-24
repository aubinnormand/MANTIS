import sys
import os
from pathlib import Path
import ipynbname
import time
import geopandas as gpd
import pandas as pd

# --- Ajouter dossier utils au Python path ---
base_path = Path(ipynbname.path()).parent.parent.parent  # racine du projet
utils_path = base_path / 'Code' / 'Utils'
if str(utils_path.resolve()) not in sys.path:
    sys.path.append(str(utils_path.resolve()))

from utils import generate_taxo_dict

def assign_grid_to_points(df_country, grid, cle_geo, sample_ratio=0.01):
    """
    Assigne à chaque point de df_country la cellule de la grille correspondante.
    """
    df_new = df_country.copy()
    
    gdf_points = gpd.GeoDataFrame(
        df_new,
        geometry=gpd.points_from_xy(df_new['lon'], df_new['lat']),
        crs=grid.crs
    )
    
    # Estimation du temps
    total_points = len(gdf_points)
    sample_size = max(1, int(total_points * sample_ratio))
    sample_points = gdf_points.sample(sample_size, random_state=42)

    start_time = time.time()
    gpd.sjoin(sample_points, grid[[cle_geo,'geometry']], how='left', predicate='within')
    elapsed_sample = time.time() - start_time
    estimated_total_time = elapsed_sample / sample_ratio
    print(f"⏱ Temps estimé pour l'ensemble des points : {estimated_total_time:.2f} s (basé sur {sample_size} points)")

    # Jointure spatiale complète
    start_time_full = time.time()
    gdf_joined = gpd.sjoin(gdf_points, grid[[cle_geo,'geometry']], how='left', predicate='within')
    elapsed_full = time.time() - start_time_full
    print(f"✅ Temps réel pour l'opération complète : {elapsed_full:.2f} s")

    df_new[cle_geo] = gdf_joined[cle_geo].values

    # Statistiques
    multiple_matches = gdf_joined.duplicated(subset=['geometry'], keep=False).sum()
    no_match = gdf_joined[cle_geo].isna().sum()
    print(f'➡️ {multiple_matches} points avec plusieurs cellules correspondantes')
    print(f'➡️ {no_match} points sans cellule correspondante')

    df_new = df_new.drop(columns='geometry', errors='ignore')
    return df_new

def aggregate_by_grid_species(df_country, cle_geo, dico_taxo, cle_ID='speciesID'):
    """
    Agrège les données par cellule, espèce, mois et année,
    puis fusionne avec le dictionnaire taxonomique.
    """
    colonnes_groupe = [cle_geo, cle_ID, 'month', 'year']
    
    df_agg = (
        df_country.groupby(colonnes_groupe, as_index=False)
                  .agg(
                      nombreObs=('nombreObs','sum'),
                      nombreObsUnique=('speciesID','nunique')
                  )
    )

    # Fusion avec dictionnaire taxonomique
    df_final = pd.merge(df_agg, dico_taxo, left_on=cle_ID, right_on=cle_ID, how='left')

    doublons = df_final.duplicated(subset=colonnes_groupe, keep=False).sum()
    if doublons > 0:
        print(f"⚠️ Attention : {doublons} doublons après fusion avec le dictionnaire taxonomique")

    return df_final

# -----------------------------
# Sauvegarde
# -----------------------------
def save_processed_biodiv(df_final, source, zone, cle_geo,path_data):
    output_path = path_data / source / 'processed' / zone/f"{source}_{zone}_{cle_geo}.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_final.to_csv(output_path, index=False)
    print(f"🎉 Données nettoyées sauvegardées dans : {output_path}")
    return output_path

    import pandas as pd

