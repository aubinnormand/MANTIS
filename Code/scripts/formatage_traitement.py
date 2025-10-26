import os
import pandas as pd
import geopandas as gpd
from formatage_geo import charger_donnees_geo, creer_grille_pays, generate_grid_path
from formatage_biodiv import traiter_chunks, fusionner_fichiers_par_domaine
from formatage_fusion import fusion_cells_by_obs, apply_fusion_to_biodiv
from fonctions_annexes_biodiv import *

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
# Fonction principale
# -------------------------
def traiter_pays_et_maille(
        country_name,
        grid_size_km,
        bornes_temporelles,
        path_data=r"C:\path\to\data",
        path_sig=r"C:\path\to\sig",
        cle_geo=CLE_GEO_DEFAULT,
        cle_ID="speciesKey",
        source="GBIF",
        fusion=False,
        var_obs="nombreObs",
        seuil_fusion=1000,
        methode_fusion="barycentre",
        max_distance=100):
    """
    Exécute le traitement complet pour un pays et une taille de maille donnés.
    Les fonctions nécessaires doivent être importées depuis les modules utils_*.
    """
    # 🔹 Charger les données géographiques
    world_terrestre, world_marin = charger_donnees_geo(path_sig, source=source)
    
    # 🔹 Créer les grilles
    country_grid_terrestre, country_grid_marin, country_grid_combined = creer_grille_pays(
        path_data, country_name, grid_size_km, world_terrestre, world_marin,
        cle_geo=cle_geo, source=source
    )

    # 🔹 Chemin du fichier raw GBIF
    path_fichier = os.path.join(
        path_data, source, "raw", f"{source}_{country_name}", f"extract{source}_{country_name}.csv"
    )
    
    colonnes_a_importer = [
        'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
        'verbatimScientificName', 'taxonRank', 'countryCode', 'occurrenceStatus',
        'individualCount', 'decimalLongitude', 'decimalLatitude', 'eventDate',
        'speciesKey', 'occurrenceID', 'year'
    ]
    
    # 🔹 Traitement des chunks
    chunk_number = traiter_chunks(
        path_fichier, colonnes_a_importer,
        country_grid_terrestre, country_grid_marin, country_grid_combined,
        bornes_temporelles
    )

    # 🔹 Création du dossier processed
    processed_path = os.path.join(path_data, source, "processed", f"{source}_{country_name}")
    os.makedirs(processed_path, exist_ok=True)

    # 🔹 Fusionner les fichiers par domaine
    for domaine in ['terrestre', 'marin', 'combined']:
        df_final = fusionner_fichiers_par_domaine(
            country_name, chunk_number, domaine,
            cle_geo=cle_geo, cle_ID=cle_ID,
            bornes_temporelles=bornes_temporelles,
            path_data=path_data, source=source
        )
    
        # --- Fusion si demandée ---
        if fusion and domaine == 'terrestre':
            grid = country_grid_terrestre
            df_obs = df_final.groupby(cle_geo, as_index=False)[var_obs].sum()
            
            merged = grid.merge(
                df_obs,
                on=cle_geo,
                how="left"
            )

            merged[var_obs] = merged[var_obs].fillna(0).astype(float)

            print(f"⚡ Application de la fusion des cellules ({methode_fusion})...")
            merged_fusion = fusion_cells_by_obs(
                merged,
                cle_geo=cle_geo, var_obs=var_obs,
                seuil=seuil_fusion,
                methode=methode_fusion,
                max_distance=max_distance
            )

            # Sauvegarde du GeoJSON fusionné
            grid_fusion_path = generate_grid_path(
                country_name, grid_size_km, domaine, cle_geo,
                path_data=path_data, source=source,
                fusion=True, seuil_fusion=seuil_fusion
            )
            merged_fusion.to_file(grid_fusion_path, driver="GeoJSON")
            print(f"✅ grid_fusion sauvegardé dans : {grid_fusion_path}")
            
            # Appliquer la fusion aux données de biodiversité
            df_maille_espece = apply_fusion_to_biodiv(
                df_final, merged_fusion,
                cle_geo=cle_geo, cle_ID=cle_ID, var_obs=var_obs
            )
            
            final_file_path = os.path.join(
                path_data, source, 'processed', f"{source}_{country_name}",
                f"data_{source}_{country_name}_{domaine}_{cle_geo}Adapt{seuil_fusion}_{cle_ID}_periodes{'_'.join(map(str,bornes_temporelles))}.csv"
            )
            df_maille_espece.to_csv(final_file_path, index=False)

    print(f"✅ Traitement complet pour {country_name} avec une taille de maille {grid_size_km} km")
    return df_final
