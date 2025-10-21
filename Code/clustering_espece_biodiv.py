#clustering_espece_biodiv.py
# Cluster par espèces

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster
from fonctions_annexes_biodiv import generer_dictionnaire_taxonomie
import matplotlib.pyplot as plt
import pandas as pd
from normalisation_biodiv import normaliser_log
import time
from scipy.spatial.distance import squareform
import numpy as np
import skfuzzy as fuzz



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

# Crée le tableau pivot mailles × espèces
def creer_pivot(df, cle_geo='codeMaille10Km', cle_ID='cdRef', col_val='nombreObs'):
    return df.pivot_table(index=cle_geo, columns=cle_ID, values=col_val, fill_value=0)


def clustering_hierarchique(df_fil, methode='ward', level=5, crit='maxclust',
                            cle_geo='codeMaille10Km', cle_ID='cdRef', col_val='nombreObs',
                            corr_method='pearson'):
    # Tableau maille × espèces
    pivot_df = creer_pivot(df_fil, cle_geo, cle_ID, col_val)

    # Matrice de corrélation
    corr_matrix = pivot_df.T.corr(method=corr_method)
    corr_matrix = corr_matrix.fillna(0)

    # Distance = 1 - corrélation
    dist_matrix = 1 - corr_matrix

    # Forcer la diagonale à 0
    np.fill_diagonal(dist_matrix.values, 0)

    # Conversion en forme condensée
    dist_condensed = squareform(dist_matrix.values)

    # Clustering hiérarchique
    Z = linkage(dist_condensed, method=methode)
    clusters = fcluster(Z, t=level, criterion=crit)

    return pd.DataFrame({
        cle_ID: corr_matrix.index,
        f'Cluster_{methode.capitalize()}': clusters
    })

def clustering_fuzzy(df_fil, n_clusters=5, m=2,
                     cle_geo='codeMaille10Km', cle_ID='cdRef', col_val='nombreObs'):
    """
    Clustering flou (Fuzzy C-Means) pour des communautés d'espèces.
    Retourne un DataFrame avec l'espèce et son cluster fuzzy.
    """

    # Créer le pivot : lignes = mailles, colonnes = espèces
    pivot_df = creer_pivot(df_fil, cle_geo, cle_ID, col_val)  # maille x espèces

    # Transposer : lignes = espèces, colonnes = mailles
    X = pivot_df.T.values

    # Supprimer les espèces entièrement nulles (jamais observées)
    non_empty = X.sum(axis=1) > 0
    X = X[non_empty, :]
    cols_kept = pivot_df.columns[non_empty]

    if X.shape[0] == 0:
        raise ValueError("Aucune espèce observée après filtrage des colonnes vides.")

    # Fuzzy C-Means
    u, _, _, _, _, _, _ = fuzz.cluster.cmeans(
        X, c=n_clusters, m=m, error=0.005, maxiter=1000, init=None
    )

    # Attribution du cluster dominant à chaque espèce
    labels = np.argmax(u, axis=0) + 1  # +1 pour commencer à 1

    # Vérification que la longueur correspond
    if len(labels) != len(cols_kept):
        raise ValueError("Longueur des labels différente du nombre d'espèces conservées.")

    # Créer le DataFrame final
    df_clusters = pd.DataFrame({
        cle_ID: cols_kept.astype(str),  # forcer en str pour merge
        'Cluster_Fuzzy': labels
    })

    return df_clusters



def clustering_gmm(df_fil, n_clusters=5,
                   cle_geo='codeMaille10Km', cle_ID='cdRef', col_val='nombreObs'):
    from sklearn.mixture import GaussianMixture
    from sklearn.preprocessing import StandardScaler
    pivot_df = creer_pivot(df_fil, cle_geo, cle_ID, col_val)
    X = pivot_df.T.values
    X_scaled = StandardScaler().fit_transform(X)
    gmm = GaussianMixture(n_components=n_clusters, random_state=42)
    gmm.fit(X_scaled)
    labels = gmm.predict(X_scaled) + 1
    return pd.DataFrame({cle_ID: pivot_df.columns,
                         'Cluster_GMM': labels})


def clustering_reseau(df_fil, seuil=0.3,
                      cle_geo='codeMaille10Km', cle_ID='cdRef', col_val='nombreObs'):
    import networkx as nx
    import community as community_louvain
    pivot_df = creer_pivot(df_fil, cle_geo, cle_ID, col_val)
    corr_matrix = pivot_df.T.corr()
    G = nx.Graph()
    for sp1 in corr_matrix.columns:
        for sp2 in corr_matrix.columns:
            if sp1 != sp2 and corr_matrix.loc[sp1, sp2] > seuil:
                G.add_edge(sp1, sp2, weight=corr_matrix.loc[sp1, sp2])
    partition = community_louvain.best_partition(G, weight='weight')
    df_louvain = pd.DataFrame.from_dict(partition, orient='index',
                                        columns=['Cluster_Louvain'])
    df_louvain.index.name = cle_ID
    return df_louvain.reset_index()


def clustering_lda(df_fil, n_clusters=5,
                   cle_geo='codeMaille10Km', cle_ID='cdRef', col_val='nombreObs'):
    from sklearn.decomposition import LatentDirichletAllocation
    pivot_df = creer_pivot(df_fil, cle_geo, cle_ID, col_val)
    lda = LatentDirichletAllocation(n_components=n_clusters, random_state=42)
    lda.fit(pivot_df.values)
    species_topic = lda.components_.T
    species_topic = species_topic / species_topic.sum(axis=1, keepdims=True)
    labels = np.argmax(species_topic, axis=1) + 1
    return pd.DataFrame({cle_ID: pivot_df.columns,
                         'Cluster_LDA': labels})


def former_communautes_all(df_fil, cle_geo='codeMaille10Km', cle_ID='cdRef', col_val='nombreObs',
                           n_clusters=5, seuil_corr=0.3, m_fuzzy=2):
    
    # Forcer cle_ID en str dans df_fil
    df_fil = df_fil.copy()
    df_fil[cle_ID] = df_fil[cle_ID].astype(str)
    
    # Créer le df final
    df_final = pd.DataFrame({cle_ID: df_fil[cle_ID].unique()})
    
    temps = {}

    # --- Ward
    print("▶ Ward...")
    t0 = time.time()
    df_ward = clustering_hierarchique(df_fil, methode='ward', level=n_clusters,
                                      cle_geo=cle_geo, cle_ID=cle_ID, col_val=col_val)
    df_ward[cle_ID] = df_ward[cle_ID].astype(str)
    df_final = df_final.merge(df_ward, on=cle_ID, how='left')
    temps['Ward'] = round(time.time() - t0, 3)

    # --- Fuzzy
    print("▶ Fuzzy...")
    t0 = time.time()
    df_fuzzy = clustering_fuzzy(df_fil, n_clusters=n_clusters, m=m_fuzzy,
                                cle_geo=cle_geo, cle_ID=cle_ID, col_val=col_val)
    df_fuzzy[cle_ID] = df_fuzzy[cle_ID].astype(str)
    df_final = df_final.merge(df_fuzzy, on=cle_ID, how='left')
    temps['Fuzzy'] = round(time.time() - t0, 3)

    # --- GMM
    print("▶ GMM...")
    t0 = time.time()
    df_gmm = clustering_gmm(df_fil, n_clusters=n_clusters,
                            cle_geo=cle_geo, cle_ID=cle_ID, col_val=col_val)
    df_gmm[cle_ID] = df_gmm[cle_ID].astype(str)
    df_final = df_final.merge(df_gmm, on=cle_ID, how='left')
    temps['GMM'] = round(time.time() - t0, 3)

    # --- Louvain
    print("▶ Louvain...")
    t0 = time.time()
    df_louvain = clustering_reseau(df_fil, seuil=seuil_corr,
                                   cle_geo=cle_geo, cle_ID=cle_ID, col_val=col_val)
    df_louvain[cle_ID] = df_louvain[cle_ID].astype(str)
    df_final = df_final.merge(df_louvain, on=cle_ID, how='left')
    temps['Louvain'] = round(time.time() - t0, 3)

    # --- LDA
    print("▶ LDA...")
    t0 = time.time()
    df_lda = clustering_lda(df_fil, n_clusters=n_clusters,
                            cle_geo=cle_geo, cle_ID=cle_ID, col_val=col_val)
    df_lda[cle_ID] = df_lda[cle_ID].astype(str)
    df_final = df_final.merge(df_lda, on=cle_ID, how='left')
    temps['LDA'] = round(time.time() - t0, 3)

    print("\n⏱ Temps de calcul (s) :")
    for k, v in temps.items():
        print(f" - {k:<8}: {v:.3f} s")

    return df_final, temps


   
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


