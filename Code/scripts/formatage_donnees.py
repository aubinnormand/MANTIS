# formatage_donnees.py
import time
import pandas as pd
import re
from datetime import datetime
from fonctions_annexes_biodiv import generer_dictionnaire_taxonomie

def filtrer_donnees_brutes(df,groupe,cle_ID='cdRef'):
    start_time = time.time()
    # Lire uniquement les colonnes spécifiées du fichier dans un DataFrame
    df = df.rename(columns={'nomScientifiqueRef': 'nomScientifique'})
    
    n_especes_entrée=len(df['cdNom'].unique())
    n_obs_entrée=len(df)
    print("En entrée : nombre d'espèces observées :",n_especes_entrée)
    print("En entrée : nombre d'obs :",n_obs_entrée)
    colonnes_obligatoires=[cle_ID]
    # Supprimer les lignes ou une donnée cruciale manque
    df_cleaned = df.dropna(subset=colonnes_obligatoires)

    # Supprimer les lignes ou une le nom scientifique ne contient pas d'espace : c'est a dire souvent pas un nom d'espece complet
    df_cleaned = df_cleaned[df_cleaned['nomScientifique'].str.contains(' ')]

    #ajouter la colonne 'all'
    df_cleaned['all']='All'

    df_cleaned['codeMaille10Km'] = df_cleaned['codeMaille10Km'].str.split('|').str[0].str.strip()
    
    df_cleaned.rename(columns={'codeInseeDepartement': 'codeDepartement'}, inplace=True)
    
    df_cleaned['nomScientifique'] = df_cleaned['nomScientifique'].apply(nettoyer_nom_scientifique)
    
    if groupe=='reptiles':
        df_cleaned['classe'] = df_cleaned['classe'].fillna('Reptilia')
    
    # Ajouter des colonnes year, month, day et day_of_year
    df_cleaned=extraire_annee_mois_jour(df_cleaned)
    
    print(f"En sortie : nombre d'espèces observées :{len(df_cleaned[cle_ID].unique())} soit une perte de {100-(round(len(df_cleaned[cle_ID].unique())/n_especes_entrée*100))}%")
    print(f"En sortie : nombre d'obs :{len(df_cleaned)} soit une perte de {100-round(len(df_cleaned)/n_obs_entrée*100)} %") 
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"Temps d'exécution : {round(execution_time)} secondes\n")

    return df_cleaned
    
def formater_maille_espece(df,cle_geo='codeMaille10Km',cle_ID='cdRef',annee_min=None,bornes_temporelles=None):
    # Convertir la colonne 'year' en int
    df['year'] = pd.to_numeric(df['year'], errors='coerce').fillna(0).astype(int)

    # Convertir la colonne 'cdNom' en int
    df[cle_ID] = df[cle_ID].astype(int)

    if annee_min is not None:
        df=df[df['year']>=annee_min]

    # Choisir des bornes temporelles et assigner une période aux données
    if bornes_temporelles is not None:
        df = df.copy()
        df.loc[:, 'periode'] = pd.cut(df['year'], bins=bornes_temporelles, 
                       labels=[f'Période {i+1}: {bornes_temporelles[i]+1} à {bornes_temporelles[i+1]}' for i in range(len(bornes_temporelles) - 1)],
                       include_lowest=False)  # include_lowest=True inclut la borne inférieureure
        # Compter le nombre de données dans chaque intervalle
        compte_par_periode = df['periode'].value_counts()
        print(compte_par_periode)
         # Compter les occurrences d'observation de chaque taxon pour chaque code et période
        df_maille_espece = df.groupby([cle_geo,cle_ID,'periode'], observed=True).size().reset_index(name='nombreObs')
       
    else:
        # Compter les occurrences d'observation de chaque taxon pour chaque code
        df_maille_espece = df.groupby([cle_geo,cle_ID], observed=True).size().reset_index(name='nombreObs')
        
    df_dico=generer_dictionnaire_taxonomie(df,cle_ID)
    df_maille_espece=pd.merge(df_maille_espece,df_dico,on=cle_ID)

    print("\n")
    
    return df_maille_espece

def extraire_code_departement(s):
    numbers = re.findall(r'\d+', s)
    numbers = [num for num in numbers if int(num) <= 900]  # Filtrer les valeurs <= 900
    if '2A' in s or '2B' in s:
        if s=="2A | 2B":
            return '2A','2B'
        if '2A' in s:
            return '2A'
        if '2B' in s:
            return '2B'
    else:
        return ', '.join(numbers)

def extraire_annee_mois_jour(df_raw,cle_date='dateObservation'):
    # Convertir la colonne 'dateObservation' en type datetime
    df_raw[cle_date] = pd.to_datetime(df_raw[cle_date],format='ISO8601', errors='coerce')

    # Créer les colonnes 'year', 'month' et 'day'
    df_raw['year'] = df_raw[cle_date].dt.year
    df_raw['month'] = df_raw[cle_date].dt.month
    df_raw['day'] = df_raw[cle_date].dt.day
    df_raw['day_of_year'] = df_raw[cle_date].dt.dayofyear
    return df_raw

def nettoyer_nom_scientifique(nom):
    mots = nom.split()  # Séparer les mots par espace
    if len(mots) == 3 and mots[1].lower() == mots[2].lower():  # Vérifier si le deuxième et troisième mot sont identiques
        return ' '.join(mots[:2])  # Garder seulement les deux premiers mots
    return nom  # Sinon, renvoyer le nom tel quel



