#clustering_espece_biodiv.py
# Cluster par espèces

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from fonctions_annexes_biodiv import generer_dictionnaire_taxonomie
import matplotlib.pyplot as plt
import pandas as pd
from normalisation_biodiv import normaliser_log

def generer_dendogram(correlation_matrix,methode='ward',display=0):
    # Appliquer le clustering hiérarchique sur la matrice de corrélation
    Z = linkage(correlation_matrix, method=methode)
    
    # Visualiser le dendrogramme
    if display==1:
        plt.figure(figsize=(10, 7))
        dendrogram(Z, labels=correlation_matrix.columns)
        plt.show()
    return Z

def former_cluster_espece(df_input, Z, col_valeur='nombreObs',cle_geo='codeMaille10Km',cle_ID='cdRef', level=4, crit='distance'):
    # Créer un tableau croisé à partir du DataFrame d'entrée
    pivot_df = df_input.pivot_table(index=cle_geo, columns=cle_ID, values=col_valeur, fill_value=0)

    # Choisir un seuil pour définir les clusters
    clusters = fcluster(Z, t=level, criterion=crit)

    # Transposer le tableau croisé pour ajouter les clusters
    pivot_df_trans = pivot_df.transpose()
    pivot_df_trans['Cluster_corr'] = clusters
    pivot_df_trans = pivot_df_trans.reset_index()

    # Créer le dictionnaire d'espèces
    df_dico = generer_dictionnaire_taxonomie(df_input,cle_ID)  # Vérifiez que df_input contient toutes les colonnes nécessaires

    # Fusionner les données
    df_cluster = pd.merge(pivot_df_trans[[cle_ID,'Cluster_corr']], df_dico, on=cle_ID, how='left')  # Utilisez 'how=left' pour éviter les problèmes d'index
  
    return df_cluster

def chercher_numcluster_espece(df_corr_cluster,cle_ID,sujet):
    num_cluster=df_corr_cluster[(df_corr_cluster[cle_ID]==sujet)][['Cluster_corr']]
    num_cluster=int(num_cluster.iat[0, 0])
    print(f'Le cluster contenant {sujet} est le cluster n° {num_cluster}')
    return num_cluster

def lister_especes_dans_cluster(df_input_clustered,num_cluster,col_valeur='nombreObs',cle_geo='codeMaille10Km',cle_ID='cdRef'):

    df_input_clustered['Cluster_corr'] = df_input_clustered['Cluster_corr'].astype(str).str.strip()

    df_filt = df_input_clustered[df_input_clustered['Cluster_corr'] == str(num_cluster)]

    liste_especes=df_filt[cle_ID].unique()

    grouped = df_filt.groupby([cle_ID])[col_valeur].sum()

    dico_taxo=generer_dictionnaire_taxonomie(df_filt,cle_ID)
    dico_taxo = dico_taxo.reset_index(drop=True)
    grouped=pd.merge(grouped,dico_taxo,on=cle_ID)
    
    # Sort the individual taxon sums in descending order
    grouped = grouped.sort_values(by=col_valeur,ascending=False)
    liste_espece_cluster = grouped.reset_index(drop=True)

    return liste_espece_cluster

def lister_espece_cluster(df_input,df_corr_cluster,num_cluster,col_valeur='nombreObs'):
    
    liste_cluster=afficher_liste_cluster(df_cluster_espece,df_input,col_valeur)
    return liste_cluster

def grouper_dans_maille(df_input, df_corr_cluster, col_valeur='nombreObs', cle_geo='codeMaille10Km', 
                        cle_ID='cdRef', cle_ID_nomScientifique='species', cleID_nomVernaculaire='vernacularName_fr'):

    # Associer les clusters aux espèces
    df_input_cluster_espece = df_input.merge(df_corr_cluster[[cle_ID, 'Cluster_corr']], on=cle_ID)

    # Calculer le nombre d'espèces par cluster
    df_nombre_espece = df_corr_cluster.groupby('Cluster_corr')[cle_ID].nunique().rename('nombre_espece').reset_index()

    # Regrouper les observations par maille et cluster
    colonnes_obs = [col for col in df_input.columns if 'nombreObs' in col]
    grouped = df_input_cluster_espece.groupby([cle_geo, 'Cluster_corr'])[colonnes_obs].sum().reset_index()

    # Normalisation par le nombre d'espèces dans chaque cluster
    grouped = grouped.merge(df_nombre_espece, on='Cluster_corr')
    for col in ['nombreObs_norm_par_espece', col_valeur, 'nombreObs']:
        grouped[col] = grouped[col] / grouped['nombre_espece']

    # Ajouter les noms scientifiques et vernaculaires des espèces par cluster
    df_filt = df_input_cluster_espece.groupby([cle_ID, 'Cluster_corr'])[col_valeur].sum().reset_index()
    df_dico = generer_dictionnaire_taxonomie(df_input, cle_ID)
    df_filt = df_filt.merge(df_dico, on=cle_ID, how='left').sort_values(by=col_valeur, ascending=False)

    grouped_noms = df_filt.groupby('Cluster_corr')[[cle_ID_nomScientifique, cleID_nomVernaculaire]].agg(
        lambda x: ', '.join(x.dropna().unique())).reset_index()
    
    grouped = grouped.merge(grouped_noms, on='Cluster_corr').sort_values(by=col_valeur, ascending=False)

    # Réordonner les colonnes
    colonnes_ordre = [cle_geo, 'Cluster_corr', cle_ID_nomScientifique, cleID_nomVernaculaire] + colonnes_obs
    return grouped[colonnes_ordre].reset_index(drop=True)


def etude_cluster_local(df_input,liste_mailles,df_corr_cluster,num_cluster,col_values,cle_ID,cle_geo):

    df_filt=df_input[df_input[cle_geo].isin(liste_mailles)]
    
    liste_espece_cluster_local=lister_especes_dans_cluster(df_filt,num_cluster,col_values,cle_geo,cle_ID)
    
    liste_espece_cluster_entier=lister_especes_dans_cluster(df_input,num_cluster,col_values,cle_geo,cle_ID)

    liste_espece_merged= pd.merge(liste_espece_cluster_local, liste_espece_cluster_entier,how='outer')

    liste_espece_merged = liste_espece_merged.sort_values(by=col_values, ascending=False)
    liste_espece_merged=liste_espece_merged.reset_index(drop=True)
    return liste_espece_merged

def chercher_especes_pas_presentes(df_input,df_cluster_maille,df_corr_cluster,liste_codes,col_choice='nombreObs_normalisé_par_espece'):

    # chercher les espèces qui pourraient être potentiellement présentes
    liste_cluster=df_cluster_maille['Cluster_corr'].unique()[:20]
    df_espece_pas_presentes=pd.DataFrame()
    for num_cluster in liste_cluster:
        df_cluster_local=etude_cluster_local(df_input,liste_codes,df_corr_cluster,num_cluster,col_choice)
        df_espece_pas_presentes_temp=df_cluster_local[df_cluster_local['nombreObs']==0]
        df_espece_pas_presentes=pd.concat([df_espece_pas_presentes, df_espece_pas_presentes_temp], ignore_index=True)
    df_espece_pas_presentes=pd.merge(df_espece_pas_presentes,df_corr_cluster[['nomScientifique','Cluster_corr']],on='nomScientifique')
    colonnes_obs = [col for col in df_input.columns if 'nombreObs' in col]
    df_espece_pas_presentes=df_espece_pas_presentes.drop(columns=colonnes_obs)

    df_input_grouped=df_input.groupby(['nomScientifique'])[colonnes_obs].sum()
    df_input_grouped=df_input_grouped.reset_index(drop=False)
    df_espece_pas_presentes=pd.merge(df_espece_pas_presentes,df_input_grouped[['nomScientifique','nombreObs','nombreObs_normalisé_par_maille_regne']],on='nomScientifique')
    df_espece_pas_presentes = df_espece_pas_presentes.groupby('Cluster_corr', sort=False).apply(lambda x: x.sort_values('nombreObs', ascending=False)).reset_index(drop=True)
    return df_espece_pas_presentes

def chercher_zone_espece_pas_presente(df_input,df_corr_cluster,df_inpn_cluster_espece,espece,
                                      col_choice='nombreObs_normalisé_par_maille_regne'):

    num_cluster=chercher_cluster_espece(df_corr_cluster,df_input,col_choice,espece)
    
    df_cluster_espece_manquante=df_inpn_cluster_espece[df_inpn_cluster_espece['Cluster_corr']==num_cluster]
    df_cluster_espece_manquante = df_cluster_espece_manquante.sort_values(by=col_choice,ascending=False)
    
    return df_cluster_espece_manquante
