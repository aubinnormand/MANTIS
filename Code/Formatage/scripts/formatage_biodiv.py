# biodiv_utils.py

import pandas as pd
import os, csv
from pathlib import Path
import numpy as np
import importlib

import formatage_biodiv_config
importlib.reload(formatage_biodiv_config)
from formatage_biodiv_config import colonnes_import, col_map

# -----------------------------
# Fonction utilitaire pour afficher la réduction
# -----------------------------
def _log_reduction(df_before, df_after, step_name):
    """Affiche la perte de données entre deux étapes."""
    n_before = len(df_before)
    n_after = len(df_after)
    removed = n_before - n_after
    pct = 0 if n_before == 0 else round((removed / n_before) * 100, 2)
    print(f"   🔧 {step_name:<30} -{removed} lignes ({pct}%)")

# -----------------------------
# Détecteur de séparateur
# -----------------------------
def detect_sep(path_fichier, seps=[",",";","\t","|"," "], nrows=5):
    best_sep = None
    best_cols = 0
    for sep in seps:
        try:
            df = pd.read_csv(path_fichier, sep=sep, nrows=nrows, engine="python")
            if len(df.columns) > best_cols:
                best_cols = len(df.columns)
                best_sep = sep
        except:
            continue
    if best_sep is None:
        raise ValueError("Impossible de déterminer un séparateur correct")
    print(f"✅ Séparateur choisi : '{best_sep}' avec {best_cols} colonnes")
    return best_sep

# -----------------------------
# Nettoyage des données
# -----------------------------
def clean_biodiv_data(df, cle_ID, annee_min=1):
    df = df.copy()
    print(f"🧾 Colonnes présentes : {list(df.columns)}")

    # ---------------- Dates ----------------
    before = df.copy()
    df['eventDate'] = pd.to_datetime(df.get('eventDate', pd.NaT), errors='coerce', utc=True)
    df['year'] = pd.to_numeric(df['year'].fillna(df['eventDate'].dt.year), errors='coerce')
    df['month'] = pd.to_numeric(df['month'].fillna(df['eventDate'].dt.month), errors='coerce')
    df = df.dropna(subset=['year', 'month']).reset_index(drop=True)
    df['year'] = df['year'].astype(int)
    df['month'] = df['month'].astype(int)
    _log_reduction(before, df, "Dates (année + mois valides)")

    # ---------------- Coordonnées ----------------
    before = df.copy()
    df['lon'] = pd.to_numeric(df.get('lon', np.nan), errors='coerce')
    df['lat'] = pd.to_numeric(df.get('lat', np.nan), errors='coerce')
    df = df.dropna(subset=['lon', 'lat']).reset_index(drop=True)
    _log_reduction(before, df, "Coordonnées valides")

    # ---------------- Occurrence PRESENT ----------------
    before = df.copy()
    df = df[df['occurrenceStatus'] == 'PRESENT'].reset_index(drop=True)
    _log_reduction(before, df, "occurrenceStatus=PRESENT")

    # ---------------- Filtrer années min ----------------
    if annee_min is not None:
        before = df.copy()
        df = df[df['year'] >= annee_min].reset_index(drop=True)
        _log_reduction(before, df, f"Année >= {annee_min}")

    # ---------------- Filtrer taxonRank = SPECIES ----------------
    before = df.copy()
    df = df[df['taxonRank'].str.upper().isin(['SPECIES', 'VARIETY'])].reset_index(drop=True)
    _log_reduction(before, df, "taxonRank=SPECIES|VARIETY")

    # ---------------- ID espèces ----------------
    before = df.copy()

    ## Nettoyage des espaces
    df[cle_ID] = df[cle_ID].astype(str).str.strip()
    
    # Vérifier si toute la colonne est numérique
    if df[cle_ID].str.isdigit().all():
        # Conversion en int puis en str (optionnel)
        df[cle_ID] = df[cle_ID].astype(int).astype(str)
    else:
        # On garde la colonne telle quelle en string
        df[cle_ID] = df[cle_ID].astype(str)

    df = df.dropna(subset=[cle_ID])  # ➤ supprime les lignes sans ID valide
    
    _log_reduction(before, df, "species ID valides")

    # ---------------- Individual count ----------------
    df['nombreObs'] = pd.to_numeric(df['nombreObs'], errors='coerce').fillna(1)
    #df['nombreObs'] = 1 #Pour ne pas prendre en compte individualCounts

    # ---------------- Supprimer colonnes inutiles ----------------
    for col in ['eventDate', 'occurrenceStatus', 'taxonRank']:
        if col in df.columns:
            df = df.drop(columns=col)

    return df

# -----------------------------
# Harmonisation colonnes source → standard
# -----------------------------
def harmonize_columns(df, source):
    mapping = col_map.get(source)
    if mapping:
        df = df.rename(columns=mapping)
    return df

# -----------------------------
# Lecture par chunks
# -----------------------------
def read_biodiv_chunks(
    path_fichier,
    source,
    cle_ID,
    zone,
    path_data,
    annee_min=1950,
    chunksize=10_000_000
):
    colonnes_a_importer = colonnes_import.get(source)
    if colonnes_a_importer is None:
        raise ValueError(f"Source {source} non reconnue")

    sep = detect_sep(path_fichier)
    chunk_paths = []
    chunk_number = 0

    for df_chunk in pd.read_csv(
        path_fichier,
        sep=sep,
        usecols=lambda c: c in colonnes_a_importer,
        chunksize=chunksize,
        on_bad_lines='skip'
    ):
        chunk_number += 1
        print(f"\n--- Traitement chunk {chunk_number} ---")

        df_chunk = harmonize_columns(df_chunk, source)
        df_clean = clean_biodiv_data_source(
            df_chunk,
            source=source,
            annee_min=annee_min
        )

        chunk_path = save_clean_chunk(
            df_clean,
            source=source,
            zone=zone,
            path_data=path_data,
            chunk_number=chunk_number
        )

        chunk_paths.append(chunk_path)

    print("\n🔗 Concaténation des chunks...")
    df_final = pd.concat(
        (pd.read_csv(p) for p in chunk_paths),
        ignore_index=True
    )

    final_path = save_clean_biodiv(df_final, source, zone, path_data)

    print("🧹 Suppression des fichiers intermédiaires...")
    for p in chunk_paths:
        p.unlink()

    print(f"🎉 Pipeline terminé : {len(df_final)} lignes")
    return df_final

# -----------------------------
# Sauvegarde
# -----------------------------
def save_clean_chunk(df_chunk, source, zone, path_data, chunk_number):
    """
    Sauvegarde un chunk nettoyé avec suffixe _chunkXX
    """
    output_dir = Path(path_data) / 'Biodiv' / source / 'clean' / 'chunks'
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{source}_{zone}_chunk{chunk_number:03d}.csv"
    df_chunk.to_csv(output_path, index=False)

    print(f"💾 Chunk {chunk_number} sauvegardé : {output_path.name}")
    return output_path


def save_clean_biodiv(df_biodiv_clean, source, zone, path_data):
    output_path = Path(path_data) / 'Biodiv'/source / 'clean' / f"{source}_{zone}.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_biodiv_clean.to_csv(output_path, index=False)
    print(f"🎉 Données nettoyées sauvegardées dans : {output_path}")
    return output_path

    import pandas as pd

def inspect_columns_with_examples(path_fichier):
    """Affiche les colonnes et un exemple de valeur pour chacune."""
    
    sep = detect_sep(path_fichier)  # utilise ton détecteur existant
    df_tmp = pd.read_csv(path_fichier, sep=sep, nrows=200, on_bad_lines="skip")  # petite lecture
    
    print(f"\n🔎 Aperçu du fichier : {path_fichier.name}")
    print(f"📦 Colonnes détectées : {len(df_tmp.columns)}\n")
    
    for col in df_tmp.columns:
        example = df_tmp[col].dropna().iloc[0] if df_tmp[col].dropna().size > 0 else "❌ NA uniquement"
        print(f"   • {col:<25} → Exemple : {example}")
    
    return df_tmp

def _clean_dates(df, date_col='eventDate', year_col='year', month_col='month'):
    before = df.copy()
    df[date_col] = pd.to_datetime(df.get(date_col, pd.NaT), errors='coerce', utc=True)
    df[year_col] = pd.to_numeric(df.get(year_col, df[date_col].dt.year), errors='coerce')
    df[month_col] = pd.to_numeric(df.get(month_col, df[date_col].dt.month), errors='coerce')
    df = df.dropna(subset=[year_col, month_col]).reset_index(drop=True)
    df[year_col] = df[year_col].astype(int)
    df[month_col] = df[month_col].astype(int)
    _log_reduction(before, df, "Dates (année + mois valides)")
    return df

def _clean_coords(df, lon='lon', lat='lat'):
    before = df.copy()
    df[lon] = pd.to_numeric(df.get(lon, np.nan), errors='coerce')
    df[lat] = pd.to_numeric(df.get(lat, np.nan), errors='coerce')
    df=df.dropna(subset=[lon, lat]).reset_index(drop=True)
    _log_reduction(before, df, "clean coords")
    return df

def _clean_occurrence(df):
    before = df.copy()
    df=df[df['occurrenceStatus'] == 'PRESENT'].reset_index(drop=True)
    _log_reduction(before, df, "Occurencestatuts")
    return df

def _clean_ID(df, cle_ID):
    before = df.copy()
    # Nettoyage robuste des IDs
    df[cle_ID] = (
        df[cle_ID]
        .astype(str)
        .str.strip()
        .replace({"nan": None, "None": None})
    )
    # Conversion numérique si possible (gère les .0)
    df[cle_ID] = pd.to_numeric(df[cle_ID], errors="coerce")
    # Retour en string propre
    df[cle_ID] = df[cle_ID].apply(
        lambda x: str(int(x)) if pd.notna(x) else x
    )
    # Suppression des NaN
    df = df.dropna(subset=[cle_ID])
    _log_reduction(before, df, "clean ID")
    return df


def _clean_nombreObs(df):
    before = df.copy()
    #df['nombreObs'] = pd.to_numeric(df.get('nombreObs', 1), errors='coerce').fillna(1)
    df['nombreObs'] = 1
    _log_reduction(before, df, "clean nombreObs")
    return df
    
def _clean_anneemin(df,annee_min=1950):
    before = df.copy()
    df = df[df['year'] >= annee_min].reset_index(drop=True)
    _log_reduction(before, df, f"Année >= {annee_min}")
    return df

def _clean_taxonRank(df,var='taxonRank',rank=['SPECIES', 'VARIETY']):
    before = df.copy()
    #df['nombreObs'] = pd.to_numeric(df.get('nombreObs', 1), errors='coerce').fillna(1)
    df=df[df[var].str.upper().isin(rank)].reset_index(drop=True)
    _log_reduction(before, df, "clean taxonRank")
    return df

def clean_gbif(df,annee_min=1950):
    df = df.copy()
    df = _clean_dates(df, 'eventDate')
    df = _clean_anneemin(df,annee_min)
    df = _clean_coords(df, lon='lon', lat='lat')
    df = _clean_occurrence(df)
    df = _clean_taxonRank(df,var='taxonRank',rank=['SPECIES', 'VARIETY'])
    df = _clean_ID(df, 'speciesID')
    df = _clean_nombreObs(df)
    # supprimer colonnes inutiles
    for col in ['eventDate','occurrenceStatus','taxonRank']:
        if col in df.columns: df = df.drop(columns=col)
    return df
    
        
def clean_inpn(df,annee_min=1950):
    df = df.copy()
    df = _clean_dates(df, 'eventDate')
    df = _clean_anneemin(df,annee_min)
    df = _clean_coords(df, lon='lon', lat='lat')
    df = _clean_occurrence(df)
    df = df[df['taxonRank'].str.upper().isin(['SPECIES', 'VARIETY'])].reset_index(drop=True)
    df = _clean_ID(df, 'cdRef')
    df = _clean_nombreObs(df)
    for col in ['dateObs','occurrenceStatus','taxonRank']:
        if col in df.columns: df = df.drop(columns=col)
    return df

def clean_silene(df,annee_min=1950):
    df = df.copy()
    df = _clean_dates(df, 'eventDate')
    df = _clean_anneemin(df,annee_min)
    df = _clean_coords(df, lon='lon', lat='lat')
    # species = deux premiers mots
    if 'nom_valide' in df.columns:
        df['species'] = df['nom_valide'].astype(str).str.split().str[:2].str.join(' ')
        df['genus'] = df['nom_valide'].astype(str).str.split().str[0]
    df = _clean_ID(df, 'speciesID')
    df = _clean_nombreObs(df)
    
    required_cols = ['speciesID','nombreObs','species','genus','family','order','class','kingdom','month','year','lon','lat']
    missing = [c for c in required_cols if c not in df.columns]
    
    if missing:
        raise KeyError(f"❌ Colonnes manquantes dans SILENE : {missing}\n💡 Colonnes présentes : {list(df.columns)}")
    
    df = df[required_cols]

    return df

def clean_biodiv_data_source(df, source,annee_min=1950):
    if source == 'GBIF':
        return clean_gbif(df,annee_min)
    elif source == 'INPN':
        return clean_inpn(df,annee_min)
    elif source == 'SILENE':
        return clean_silene(df,annee_min)
    else:
        raise ValueError(f"Source inconnue : {source}")


