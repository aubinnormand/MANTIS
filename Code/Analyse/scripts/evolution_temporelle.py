# evolution_temporelle.py
import pandas as pd
import numpy as np
import geopandas as gpd

def suivre_disparition_geo(df,cle_ID='cdRef',cle_geo='codeMaille10Km'):
    # Filtrer pour les espèces distinctes dans chaque maille pour chaque période
    periode_1_df = df[df['periode'].str.contains('Période 1')][[cle_geo, cle_ID]].drop_duplicates()
    periode_2_df = df[df['periode'].str.contains('Période 2')][[cle_geo, cle_ID]].drop_duplicates()

    count_total=df.groupby(cle_geo).size().reset_index(name='nb_total')
    
    # Nombre d'espèces présentes dans la période 1 par maille
    count_p1 = periode_1_df.groupby(cle_geo).size().reset_index(name='nb_p1')
    
    # Nombre d'espèces présentes dans la période 2 par maille
    count_p2 = periode_2_df.groupby(cle_geo).size().reset_index(name='nb_p2')
    
    # Espèces présentes dans les deux périodes (jointure interne)
    both_periods = pd.merge(periode_1_df, periode_2_df, on=[cle_geo, cle_ID])
    count_both = both_periods.groupby(cle_geo).size().reset_index(name='nb_p1_and_p2')
    
    # Fusion des DataFrames pour obtenir toutes les informations par maille
    suivi_df = pd.merge(count_total, count_p1, on=cle_geo, how='outer').fillna(0)
    suivi_df = pd.merge(suivi_df, count_both, on=cle_geo, how='outer').fillna(0)
    suivi_df = pd.merge(suivi_df, count_p2, on=cle_geo, how='outer').fillna(0)
    
    # Calcul des espèces présentes dans période 1 mais pas période 2
    suivi_df['nb_p1_not_p2'] = suivi_df['nb_p1'] - suivi_df['nb_p1_and_p2']
    
    # Calcul des espèces présentes dans période 2 mais pas période 1
    suivi_df['nb_p2_not_p1'] = suivi_df['nb_p2'] - suivi_df['nb_p1_and_p2']
    
    suivi_df['taux_disparue'] = np.where(suivi_df['nb_p1'] != 0, 
                                     (suivi_df['nb_p1'] - suivi_df['nb_p1_and_p2']) / suivi_df['nb_p1'] * 100, 
                                     np.nan)
    suivi_df['taux_apparue'] = (suivi_df['nb_p2'] - suivi_df['nb_p1_and_p2'])/suivi_df['nb_total'] *100 
                                         

    # Finalisation du DataFrame avec les colonnes demandées
    suivi_df = suivi_df[[cle_geo, 'nb_total',
                         'nb_p1', 'nb_p1_and_p2', 
                         'nb_p1_not_p2','taux_disparue',
                         'nb_p2', 'nb_p2_not_p1','taux_apparue']]
    
    return suivi_df

def determiner_statut(periodes):
    if 'Période 1: 1801 à 1990' in periodes:
        if 'Période 2: 1991 à 2010' in periodes and 'Période 3: 2011 à 2024' in periodes:
            return "Présence continue jusqu'à aujourd'hui"
        elif 'Période 2: 1991 à 2010' not in periodes and 'Période 3: 2011 à 2024' in periodes:
            return "Colonisation 2011-2024"
        elif 'Période 2: 1991 à 2010' not in periodes and 'Période 3: 2011 à 2024' not in periodes:
            return "Disparition 1801-1990"
        elif 'Période 2: 1991 à 2010' in periodes and 'Période 3: 2011 à 2024' not in periodes:
            return "Disparition 1991-2011"

    else:
        if 'Période 2: 1991 à 2010' in periodes and 'Période 3: 2011 à 2024' in periodes:
            return "Colonisation 1991-2010"
        elif 'Période 2: 1991 à 2010' not in periodes and 'Période 3: 2011 à 2024' in periodes:
            return "Colonisation 2011-2024"
        elif 'Période 2: 1991 à 2010' in periodes and 'Période 3: 2011 à 2024' not in periodes:
            return "Présence 1991-2010 mais absent aujourd'hui"

            

    
"""
def formater_periode
df_filt=df_papillons.copy()
cle_periode='periode'

df_filt.loc[df_filt[cle_periode]=='Période 1: 1801 à 1990',cle_periode]='Période 1: 1801 à 2010'
df_filt.loc[df_filt[cle_periode]=='Période 2: 1991 à 2010',cle_periode]='Période 1: 1801 à 2010'
df_filt.loc[df_filt[cle_periode]=='Période 3: 2011 à 2024',cle_periode]='Période 2: 2011 à 2024'

#generer le dicotionnaire d'espece
dico_taxo=generer_dictionnaire_taxonomie(df_filt,cle_ID)

# permet de garder seulement les variables que l'on souhaite
df_filt = df_filt.groupby([cle_geo,cle_ID,cle_periode],as_index=False).agg({'nombreObs': 'sum'})
df_filt=pd.merge(df_filt,dico_taxo,on=cle_ID)
# Normalisations des données
df_filt=normaliser_unique(df_filt)
df_filt=normaliser_par_espece(df_filt,'cdRef')
df_filt=normaliser_par_maille_et_clade(df_filt, code_col=cle_geo, clade_col='regne', observation_col='nombreObs')
df_filt=normaliser_log(df_filt,'nombreObs_norm_par_maille_et_regne')
""
df_filt=normaliser_par_periode(df_filt)
df_filt=normaliser_par_periode(df_filt,'nombreObs_norm_par_espece')
df_filt=normaliser_par_periode(df_filt,'nombreObs_norm_par_maille_et_regne')
"""