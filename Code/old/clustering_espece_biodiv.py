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

def lister_especes_dans_cluster(df_input_clustered,num_cluster,colonnes_obs='nombreObs',col_valeur='nombreObs',cle_geo='codeMaille10Km',cle_ID='cdRef'):

    df_filt = df_input_clustered[df_input_clustered['Cluster_corr'] == num_cluster]

    liste_especes=df_filt[cle_ID].unique()

    grouped = df_filt.groupby([cle_ID])[colonnes_obs].sum()

    dico_taxo=generer_dictionnaire_taxonomie(df_filt,cle_ID)
    dico_taxo = dico_taxo.reset_index(drop=True)
    grouped=pd.merge(grouped,dico_taxo,on=cle_ID)
    
    # Sort the individual taxon sums in descending order
    grouped = grouped.sort_values(by=col_valeur,ascending=False)
    liste_espece_cluster = grouped.reset_index(drop=True)

    return liste_espece_cluster


def grouper_par_cluster(df_input, df_corr_cluster, colonnes_obs='nombreObs', col_valeur='nombreObs', cle_geo='codeMaille10Km', 
                        cle_ID='cdRef', cle_ID_nomScientifique='species', cleID_nomVernaculaire='vernacularName_fr',crit='sum',norm=True):

    # Associer les clusters aux espèces
    df_input_cluster_espece = df_input.merge(df_corr_cluster[[cle_ID, 'Cluster_corr']], on=cle_ID)

    # Calculer le nombre d'espèces par cluster
    df_nombre_espece = df_corr_cluster.groupby('Cluster_corr')[cle_ID].nunique().rename('nombre_espece').reset_index()

    # Regrouper les observations par maille et cluster
    colonnes_obs = [col for col in df_input.columns if 'nombreObs' in col]
    grouped = df_input_cluster_espece.groupby([cle_geo, 'Cluster_corr'])[colonnes_obs].sum().reset_index()

    # Normalisation par le nombre d'espèces dans chaque cluster
    if norm:
        grouped = grouped.merge(df_nombre_espece, on='Cluster_corr')
        for col in colonnes_obs:
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
    df_grouper_par_maille_local= grouped[colonnes_ordre].reset_index(drop=True)


    df_grouper_par_cluster_local = df_grouper_par_maille_local.groupby('Cluster_corr').agg(
        {
            cle_ID_nomScientifique: lambda x: ', '.join(x.unique()),
            cleID_nomVernaculaire: lambda x: ', '.join(x.unique()),
            **{col: crit for col in colonnes_obs}
        }
    ).reset_index()
    
    df_grouper_par_cluster_local = df_grouper_par_cluster_local.sort_values(by=col_valeur, ascending=False)

    return df_grouper_par_cluster_local

def etudier_un_cluster_local(df_global_cluster, df_local_cluster, cle_ID, num_cluster=1,colonnes_obs='nombreObs_norm_par_maille_et_kingdom',
                            col_choice_1='nombreObs',col_choice_2='nombreObs_norm_par_espece'):

    liste_especes_cluster_choisi = lister_especes_dans_cluster(df_global_cluster, num_cluster, colonnes_obs,col_choice_1, cle_ID=cle_ID)
    
    # Liste des espèces du cluster sélectionné, triée selon col_valeur
    liste_especes_cluster_choisi_local = lister_especes_dans_cluster(df_local_cluster, num_cluster,colonnes_obs, col_valeur=col_choice_2, cle_ID=cle_ID)
    
    nouvelle_partie = liste_especes_cluster_choisi[
        ~liste_especes_cluster_choisi[cle_ID].isin(liste_especes_cluster_choisi_local[cle_ID])
    ].copy()
    
    # Filtrer les lignes de liste_especes_cluster_choisi qui ne sont pas dans liste_especes_cluster_choisi_local
    nouvelle_partie = liste_especes_cluster_choisi[
        ~liste_especes_cluster_choisi[cle_ID].isin(liste_especes_cluster_choisi_local[cle_ID])
    ].copy()
    
    nouvelle_partie[colonnes_obs] = 0
    # Concaténer les deux DataFrames
    resultat = pd.concat([liste_especes_cluster_choisi_local, nouvelle_partie], ignore_index=True)
        
    return resultat

def rechercher_especes_localement_absentes(df_global_cluster, df_local_cluster,df_grouper_par_cluster_local, cle_ID,colonnes_obs, col_choice='nombreObs',seuil=None,max_cluster=10):
    resultat_final = []

    for num_cluster in df_grouper_par_cluster_local['Cluster_corr'].unique()[:max_cluster]:  # De 0 à 100 inclus
        resultat = etudier_un_cluster_local(df_global_cluster, df_local_cluster, cle_ID, num_cluster,colonnes_obs,'nombreObs',col_choice)
        if not resultat.empty:
            if seuil is not None:
                resultat = resultat[resultat[col_choice] < seuil]
            else:
                resultat = resultat[resultat['nombreObs'] == 0]  # Filtrer où nombreObs > 0
                
            resultat['num_cluster'] = num_cluster  # Ajouter la colonne du numéro de cluster
            resultat_final.append(resultat)

    return pd.concat(resultat_final, ignore_index=False)


