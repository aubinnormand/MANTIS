#utils.py

from pathlib import Path
import sys
import ipynbname

# ========================================================
# 🗂 Fonctions utilitaires
# ========================================================

def get_base_path(levels_up: int = 3) -> Path:
    """Retourne le chemin racine du projet depuis le notebook."""
    return Path(ipynbname.path()).parent.parents[levels_up - 1]

def print_paths(**paths):
    """Affiche proprement les chemins configurés."""
    print("📌 **CONFIGURATION DES CHEMINS**")
    for name, path in paths.items():
        print(f"📍 {name:<12}: {path}")

def ensure_folder(path: Path):
    """Crée le dossier s'il n'existe pas."""
    path.mkdir(parents=True, exist_ok=True)
    print(f"📁 Dossier créé ou existant : {path}")

def add_to_sys_path(path: Path):
    """Ajoute un dossier au Python path s'il n'est pas déjà présent."""
    resolved = str(path.resolve())
    if resolved not in sys.path:
        sys.path.append(resolved)
        print(f"✅ Ajout au sys.path : {resolved}")

# ========================================================
# Génération du dictionnaire taxonomique
# ========================================================

def generate_taxo_dict(df_country, cle_ID='speciesID'):
    """
    Génère un dictionnaire taxonomique à partir du DataFrame.
    """
    colonnes_taxo = ['speciesKey','cdRef','cdNom','taxonID','speciesID',
                    'species','vernacularName_fr','vernacularName_en','genus','family','order','class','phylum','kingdom', 
                    'nomScientifique','nomVernaculaire','genre', 'famille', 'ordre', 'classe', 'regne',
                    'especeProtegee','statutBiogeoEspeceTaxref', 'habitatEspeceTaxref','occurrenceID']
    
    colonnes_presentes = [c for c in colonnes_taxo if c in df_country.columns]

    dico_taxo = df_country[colonnes_presentes].drop_duplicates(subset=[cle_ID])
    return dico_taxo
