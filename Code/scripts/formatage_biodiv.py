import os
import geopandas as gpd
import pandas as pd
import numpy as np
import csv

# -------------------------
# Paramètres globaux
# -------------------------
CLE_ID_DEFAULT = 'speciesKey'
CLE_GEO_DEFAULT = 'codeMaille100Km'
PATH_DATA = r"C:/Users/Aubin/Documents/MANTIS/Data"  # à adapter
SOURCE = "GBIF"

def generer_dictionnaire_taxonomie(df,cle_ID='CLE_ID_DEFAULT'):
    # Définir les colonnes attendues
    colonnes_attendues = ['speciesKey','cdRef','cdNom','taxonID',
                          'species','vernacularName_fr','vernacularName_en','genus','family','order','class','phylum','kingdom', 
                          'nomScientifique','nomVernaculaire','genre', 'famille', 'ordre', 'classe', 'regne',
                      'especeProtegee','statutBiogeoEspeceTaxref', 'habitatEspeceTaxref','occurrenceID']
    
    # Vérifier quelles colonnes sont présentes dans le DataFrame
    colonnes_presentes = [col for col in colonnes_attendues if col in df.columns]
    
    df_filtered=df[colonnes_presentes]

    #Trier le DataFrame par le nombre de valeurs non nulles dans chaque ligne, en ordre décroissant
    df_sorted = df_filtered.loc[df_filtered.notnull().sum(axis=1).sort_values(ascending=False).index]
    
    # Supprimer les doublons en fonction de 'cle_ID'
    dico_taxo = df_sorted.drop_duplicates(subset=[cle_ID])
    
    # Sélectionner uniquement les colonnes présentes
    dico_taxo = dico_taxo[colonnes_presentes]

    return dico_taxo
    
def importer_biodiv(path_fichier, colonnes_a_importer=None, chunksize=10_000_000):
    """
    Générateur qui lit un fichier CSV de biodiversité par chunk.

    Args:
        path_fichier (str or Path): Chemin vers le fichier CSV.
        colonnes_a_importer (list or None): Liste des colonnes à importer. Si None, toutes les colonnes sont lues.
        chunksize (int): Taille des chunks.

    Yields:
        DataFrame: chunk du CSV filtré sur les colonnes demandées.
    """
    sep = detect_sep(path_fichier)
    for df_chunk in pd.read_csv(
        path_fichier,
        sep=sep,
        quoting=csv.QUOTE_NONE,
        chunksize=chunksize,
        on_bad_lines='skip',
        usecols=(lambda c: c in colonnes_a_importer) if colonnes_a_importer is not None else None
    ):
        yield df_chunk


def apercu_biodiv(path_fichier, colonnes_a_importer=None, n_lignes=5):
    """
    Affiche un aperçu d'un fichier de biodiversité en utilisant le premier chunk.

    Args:
        path_fichier (str or Path): Chemin vers le fichier CSV.
        colonnes_a_importer (list or None): Liste des colonnes à importer. Si None, toutes les colonnes sont lues.
        n_lignes (int): Nombre de lignes à afficher pour l'aperçu.

    Returns:
        DataFrame: Premier chunk du fichier (pour inspection).
    """
    # Récupérer le premier chunk via importer_biodiv
    chunks = importer_biodiv(path_fichier, colonnes_a_importer, chunksize=10_000_000)
    premier_chunk = next(chunks)  # récupère le premier chunk

    # Affichage
    print(f"✅ Colonnes disponibles : {list(premier_chunk.columns)}")
    print(f"✅ Nombre de lignes dans le premier chunk : {len(premier_chunk)}")
    display(premier_chunk.head(n_lignes))

    return premier_chunk
    
    def apercu_colonne(df, col="occurrenceStatus"):
        """
        Affiche le nombre de lignes pour chaque valeur unique dans la colonne occurrenceStatus.
    
        Args:
            df (pd.DataFrame): DataFrame à analyser.
            col (str): Nom de la colonne, par défaut 'occurrenceStatus'.
    
        Returns:
            pd.DataFrame: Tableau récapitulatif des valeurs et leur nombre.
        """
        counts = df[col].value_counts(dropna=False).reset_index()
        counts.columns = [col, "Nombre"]
        counts["%"] = (counts["Nombre"] / len(df) * 100).round(2)
    return counts


def process_biodiv_data(df, cle_ID=CLE_ID_DEFAULT, cle_geo=CLE_GEO_DEFAULT, annee_min=1):
    """Nettoie et formate les données de biodiversité."""
    df_cleaned = df.copy()
    
    n_especes_entree = len(df_cleaned[cle_ID].unique())
    n_obs_entree = len(df_cleaned)
    
    # Conversion dates et années
    df_cleaned['eventDate'] = pd.to_datetime(df_cleaned['eventDate'], errors='coerce', utc=True)
    df_cleaned['year'] = df_cleaned['year'].fillna(df_cleaned['eventDate'].dt.year)
    
    # Coordonnées numériques
    df_cleaned['decimalLongitude'] = pd.to_numeric(df_cleaned['decimalLongitude'], errors='coerce')
    df_cleaned['decimalLatitude'] = pd.to_numeric(df_cleaned['decimalLatitude'], errors='coerce')
    
    # Suppression lignes manquantes
    df_cleaned = df_cleaned.dropna(subset=['decimalLongitude', 'decimalLatitude', cle_ID]).reset_index(drop=True)

    # Filtrer occurrences présentes
    df_cleaned = df_cleaned[df_cleaned['occurrenceStatus'] == 'PRESENT'].reset_index(drop=True)

    # Conversion types
    df_cleaned[cle_ID] = df_cleaned[cle_ID].astype(int)
    df_cleaned['year'] = pd.to_numeric(df_cleaned['year'], errors='coerce').fillna(0).astype(int)

    if annee_min is not None:
        df_cleaned = df_cleaned[df_cleaned['year'] >= annee_min]

    # Renommer colonne de grille
    if 'grid_name' in df_cleaned.columns:
        df_cleaned.rename(columns={'grid_name': cle_geo}, inplace=True)
    df_cleaned['individualCount'] = pd.to_numeric(df_cleaned['individualCount'], errors='coerce').fillna(1)

    # Statistiques
    perte_especes = 100 - round(len(df_cleaned[cle_ID].unique()) / n_especes_entree * 100)
    perte_obs = 100 - round(len(df_cleaned) / n_obs_entree * 100)

    print(f"✅ En sortie : {len(df_cleaned[cle_ID].unique())} espèces (-{perte_especes}%)")
    print(f"✅ En sortie : {len(df_cleaned)} observations (-{perte_obs}%)")

    return df_cleaned


def add_grid_to_country(df_country, grid, cle_geo=CLE_GEO_DEFAULT):
    """Associe la maille du domaine à chaque point du DataFrame basé sur les coordonnées."""
    df_new = df_country.copy()

    bounds = grid['geometry'].bounds
    min_lons, max_lons = bounds['minx'].values, bounds['maxx'].values
    min_lats, max_lats = bounds['miny'].values, bounds['maxy'].values
    grid_names = grid[cle_geo].values

    grid_names_for_points = []
    s, n = 0, 0
    print(f"➡️ Association des mailles aux données")

    for lon, lat in zip(df_country['decimalLongitude'], df_country['decimalLatitude']):
        matching_grid = np.where((min_lons <= lon) & (lon <= max_lons) &
                                 (min_lats <= lat) & (lat <= max_lats))[0]
        if matching_grid.size > 0:
            grid_names_for_points.append(grid_names[matching_grid[0]])
            if matching_grid.size > 1: s += 1
        else:
            grid_names_for_points.append(None)
            n += 1

    print(f'{s} avec plusieurs grilles correspondantes')
    print(f'{n} sans grille correspondante')
    df_new[cle_geo] = grid_names_for_points
    return df_new


def formater_maille_espece_GBIF(df, cle_geo=CLE_GEO_DEFAULT, cle_ID=CLE_ID_DEFAULT, bornes_temporelles=None):
    """Formate le nombre d'observations par maille et par espèce (optionnellement par période)."""
    df_dico = generer_dictionnaire_taxonomie(df, cle_ID)
    df['year'] = pd.to_numeric(df['year'], errors='coerce').fillna(0).astype(int)

    if bornes_temporelles is not None:
        df.loc[:, 'periode'] = pd.cut(
            df['year'],
            bins=bornes_temporelles,
            labels=[f'Période {i+1}: {bornes_temporelles[i]+1} à {bornes_temporelles[i+1]}' 
                    for i in range(len(bornes_temporelles) - 1)],
            include_lowest=False
        )
        df_maille_espece = df.groupby([cle_geo, cle_ID, 'periode'], observed=True).size().reset_index(name='nombreObs')
    else:
        df_maille_espece = df.groupby([cle_geo, cle_ID], observed=True).size().reset_index(name='nombreObs')

    df_maille_espece = pd.merge(df_maille_espece, df_dico, on=cle_ID)
    return df_maille_espece


def process_chunk(df_biodiv, grid_terrestre=None, grid_marin=None, grid_combined=None,
                  bornes_temporelles=None, chunk_number=1,
                  cle_geo=CLE_GEO_DEFAULT, cle_ID=CLE_ID_DEFAULT,
                  path_data=PATH_DATA, country_name="", source=SOURCE):
    """Prétraitement et enregistrement des données par chunk pour chaque domaine."""
    print(f"\n🔍 Traitement du chunk {chunk_number}...")
    df_biodiv = process_biodiv_data(df_biodiv, cle_ID, cle_geo)

    dfs_domaine = {}
    for domaine, grid in zip(['terrestre', 'marin', 'combined'], [grid_terrestre, grid_marin, grid_combined]):
        if grid is not None:
            df_dom = add_grid_to_country(df_biodiv, grid, cle_geo).dropna(subset=[cle_geo]).reset_index(drop=True)
            dfs_domaine[domaine] = formater_maille_espece_GBIF(df_dom, cle_geo, cle_ID, bornes_temporelles)

    # Ajouter les noms vernaculaires
    dico_noms_vernaculaires = pd.read_csv(os.path.join(os.path.dirname(path_data), r"Taxo\TAXO_GBIF\dico_noms_vernaculaires_merged.csv"))
    for domaine, df_dom in dfs_domaine.items():
        dfs_domaine[domaine] = pd.merge(df_dom, dico_noms_vernaculaires, on=cle_ID, how="left")

    # Sauvegarde
    chaine_bornes = "_".join(map(str, bornes_temporelles)) if bornes_temporelles else "all"
    for domaine, df_dom in dfs_domaine.items():
        path_save = os.path.join(path_data, source, 'processed', f"{source}_{country_name}")
        os.makedirs(path_save, exist_ok=True)
        df_dom.to_csv(os.path.join(path_save,
                                   f"data_{source}_{country_name}_{domaine}_{cle_geo}_{cle_ID}_periodes{chaine_bornes}_{chunk_number}.csv"),
                      index=False)

    # Nettoyage mémoire
    del df_biodiv
    for df_dom in dfs_domaine.values():
        del df_dom

    print(f"\n✅ Chunk {chunk_number} traité et sauvegardé avec succès ! 🎉\n")


def detect_sep(path_fichier, seps=[",",";","\t","|"," "], nrows=5):
    """Détecte automatiquement le séparateur le plus probable."""
    best_sep, best_cols = None, 0
    for sep in seps:
        try:
            df = pd.read_csv(path_fichier, sep=sep, nrows=nrows, engine="python")
            ncols = len(df.columns)
            if ncols > best_cols:
                best_cols, best_sep = ncols, sep
        except Exception:
            continue
    if best_sep is None:
        raise ValueError("Impossible de déterminer un séparateur correct")
    print(f"✅ Séparateur choisi : '{best_sep}' avec {best_cols} colonnes")
    return best_sep


# -------------------------
# Fonction principale : traiter_chunks
# -------------------------
def traiter_chunks(path_fichier, colonnes_a_importer, grid_terrestre=None, grid_marin=None, grid_combined=None,
                   bornes_temporelles=None, path_data=PATH_DATA, source=SOURCE):
    """
    Lit les données par chunks et les traite pour chaque domaine.

    Args:
        path_fichier (str or Path): Chemin du fichier CSV.
        colonnes_a_importer (list): Colonnes à importer.
        grid_terrestre, grid_marin, grid_combined: grilles pour les domaines.
        bornes_temporelles (list): bornes temporelles pour découpage.
        path_data (str or Path): chemin racine.
        source (str): source des données (ex: GBIF).

    Returns:
        int: nombre de chunks traités.
    """
    chunk_number = 0

    for df_biodiv in importer_biodiv(path_fichier, colonnes_a_importer):
        chunk_number += 1
        process_chunk(
            df_biodiv, grid_terrestre, grid_marin, grid_combined,
            bornes_temporelles, chunk_number, path_data=path_data, source=source
        )

    return chunk_number



def fusionner_fichiers_par_domaine(country_name, chunk_number, domaine,
                                   cle_geo=CLE_GEO_DEFAULT, cle_ID=CLE_ID_DEFAULT,
                                   bornes_temporelles=None, path_data=PATH_DATA, source=SOURCE):
    """Fusionne les fichiers pour un domaine donné et génère le fichier final."""
    chaine_bornes = "_".join(map(str, bornes_temporelles)) if bornes_temporelles else "all"
    df_final = pd.DataFrame()
    fichiers_trouves = False 

    for i in range(1, chunk_number + 1):
        file_path = os.path.join(path_data, source, 'processed', f"{source}_{country_name}",
                                 f"data_{source}_{country_name}_{domaine}_{cle_geo}_{cle_ID}_periodes{chaine_bornes}_{i}.csv")
        if os.path.exists(file_path):
            fichiers_trouves = True
            df_temp = pd.read_csv(file_path)
            df_final = pd.concat([df_final, df_temp], ignore_index=True)

    if not fichiers_trouves:
        print(f"⚠️ Aucun fichier trouvé pour {country_name} - {domaine} - {cle_geo}.")
        return df_final

    # Regroupement des observations
    if 'periode' in df_final.columns:
        df_maille_espece = df_final.groupby([cle_geo, cle_ID, 'periode'], observed=True)['nombreObs'].sum().reset_index()
    else:
        df_maille_espece = df_final.groupby([cle_geo, cle_ID], observed=True)['nombreObs'].sum().reset_index()

    # Dictionnaire taxonomique
    df_dico = generer_dictionnaire_taxonomie(df_final, cle_ID)
    df_maille_espece = pd.merge(df_maille_espece, df_dico, on=cle_ID)

    # Sauvegarde
    final_file_path = os.path.join(path_data, source, 'processed', f"{source}_{country_name}",
                                   f"data_{source}_{country_name}_{domaine}_{cle_geo}_{cle_ID}_periodes{chaine_bornes}.csv")
    df_maille_espece.to_csv(final_file_path, index=False)

    # Suppression fichiers chunk
    for i in range(1, chunk_number + 1):
        file_path = os.path.join(path_data, source, 'processed', f"{source}_{country_name}",
                                 f"data_{source}_{country_name}_{domaine}_{cle_geo}_{cle_ID}_periodes{chaine_bornes}_{i}.csv")
        if os.path.exists(file_path):
            os.remove(file_path)

    return df_maille_espece
