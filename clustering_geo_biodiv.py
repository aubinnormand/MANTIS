# clustering_geo_biodiv.py 
# Cluster géographiques
import pandas as pd 
import geopandas as gpd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from fonctions_annexes_biodiv import generer_dictionnaire_taxonomie
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull
from heapq import nsmallest

def analyser_composantes_principales(df,cle_geo='codeMaille10Km',col_values='nombreObs',cle_ID='cdRef',variance_threshold=0.9,
                                     n_components=None,max_components=None):
    pivot_df = df.pivot_table(index=cle_geo, columns=cle_ID, values=col_values, fill_value=0)
   
    # Supposons que pivot_df soit ton DataFrame à normaliser
    scaler = StandardScaler()
    pivot_df_scaled = scaler.fit_transform(pivot_df)
    
    if max_components is None or max_components > len(df[cle_ID].unique()):
        max_components = len(df[cle_ID].unique())

    if n_components is None:

        pca = PCA(n_components=min(max_components,len(pivot_df)))
        principal_components = pca.fit_transform(pivot_df)
    
        # Calculer la variance expliquée cumulée
        cumulative_variance = pca.explained_variance_ratio_.cumsum()

        if (cumulative_variance >= variance_threshold).any():
            # Nombre minimum de composantes nécessaires pour atteindre le seuil
            n_components = (cumulative_variance >= variance_threshold).argmax() + 1
            print(f"Nombre de composantes nécessaires pour atteindre {variance_threshold*100:.1f}% de variance expliquée: {n_components}")
            pca = PCA(n_components=n_components)
            principal_components = pca.fit_transform(pivot_df)
            cumulative_variance = pca.explained_variance_ratio_.cumsum()

        else:
            # Si le seuil n'est pas atteint, utiliser 100 composantes ou le maximum possible
            n_components = min(max_components, pivot_df_scaled.shape[1])
            print(f"Seuil de variance non atteint. Variance cumulée de {max(cumulative_variance)*100:.1f} %. Utilisation de  {max_components} composantes.")
        
    else :
        # Recalculer la PCA avec le nombre optimal de composantes
        pca = PCA(n_components=n_components)
        principal_components = pca.fit_transform(pivot_df)
        cumulative_variance = pca.explained_variance_ratio_.cumsum()
    
    # Générer dynamiquement les noms des colonnes pour chaque composante
    column_names = [f'PC{i+1}' for i in range(n_components)]
    
    # Créer un DataFrame avec les résultats de la PCA
    df_pca = pd.DataFrame(data=principal_components, columns=column_names)
    
   # Ajouter la colonne cle_geo pour identifier les observations
    code_maille_column = pivot_df.index.to_frame(index=False)[[cle_geo]]
    df_pca = pd.concat([code_maille_column.reset_index(drop=True), df_pca], axis=1)
    df_pca.set_index(cle_geo, inplace=True)
    
    # Afficher la variance expliquée par chaque composante
    explained_variance = pca.explained_variance_ratio_
    print(f'Variance expliquée par les dix premières composante: {np.round(explained_variance[:10] * 100,1)}')
    print(f'Variance cumulée totale : {max(cumulative_variance)*100:.1f} %')
    return df_pca

def former_cluster_biogeo(df_pca,cle_geo,method='kmeans',k_cluster=5,n_init=10,display=True):

    if method=='kmeans':
        # Étape 3 : Appliquer K-means avec le nombre de clusters choisi (par exemple k=3)
        kmeans = KMeans(n_clusters=k_cluster,n_init=n_init)
        clusters = kmeans.fit_predict(df_pca)
    if method=='ward':
        ward = AgglomerativeClustering(n_clusters=k_cluster, linkage='ward')
        clusters = ward.fit_predict(df_pca)
    if method=='complete_linkage':
        complete_linkage = AgglomerativeClustering(n_clusters=k_cluster, linkage='complete')
        clusters = complete_linkage.fit_predict(df_pca)
    if method=='single_linkage':
        single_linkage = AgglomerativeClustering(n_clusters=k_cluster, linkage='single')
        clusters = single_linkage.fit_predict(df_pca)
    if method=='average_linkage':
        average_linkage = AgglomerativeClustering(n_clusters=k_cluster, linkage='average')
        clusters = average_linkage.fit_predict(df_pca)
    if method=='DBSCAN':
        dbscan = DBSCAN(eps=0.5, min_samples=5,n_clusters=k_cluster)
        clusters = dbscan.fit_predict(df_pca)
    if method=='SpectralClustering':
        clustering = SpectralClustering(n_clusters=k_cluster)
        clusters = clustering.fit_predict(df_pca)
        
    # Étape 4 : Ajouter les clusters au DataFrame
    df_cluster=df_pca.copy()
    df_cluster['Cluster'] = clusters

    # Étape 1 : Compter le nombre d'occurrences de chaque cluster
    cluster_counts = df_cluster['Cluster'].value_counts()
    
    # Étape 2 : Trier les clusters par nombre d'occurrences de manière décroissante
    sorted_clusters = cluster_counts.index
    
    # Étape 3 : Créer un dictionnaire de mappage pour réorganiser les clusters
    cluster_mapping = {old_cluster: new_cluster for new_cluster, old_cluster in enumerate(sorted_clusters, start=1)}
    
    # Réorganiser les clusters dans le DataFrame
    df_cluster['Cluster'] = df_cluster['Cluster'].map(cluster_mapping)
    
    # Compter le nombre de points dans chaque cluster
    cluster_counts = df_cluster['Cluster'].value_counts()
    
    # Afficher les résultats
    if display:
        print(cluster_counts)
    df_cluster.reset_index(inplace=True)
    # Afficher le DataFrame avec les clusters
    return df_cluster[[cle_geo, 'Cluster']]

def former_cluster_biogeo_avec_critere_spatial(df_pca, carte_maille, cle_geo, n_components, k, choix_methode, lambda_penalty,n_init=1,parameter=50,display=True):
    """Exécute l'ensemble du pipeline de clustering avec contrainte spatiale."""
    
    # Charger et fusionner les données
    df_pca_geo = pd.merge(df_pca.reset_index(), carte_maille[[cle_geo, 'geometry']], on=cle_geo)
    df_cluster = df_pca_geo.copy()
    
    # Extraire les features et les coordonnées spatiales
    features = df_pca_geo[[f"PC{i}" for i in range(1, n_components + 1)]].values
    spatial_coords = np.array([geom.centroid.coords[0] for geom in df_pca_geo.geometry])
    
    # Appliquer le clustering
    df_cluster['Cluster'] = kmeans_with_spatial_constraint(features, spatial_coords, k, choix_methode, lambda_penalty, n_init, parameter)
    
    # Réorganiser les clusters
    cluster_counts = df_cluster['Cluster'].value_counts()
    sorted_clusters = cluster_counts.index
    cluster_mapping = {old_cluster: new_cluster for new_cluster, old_cluster in enumerate(sorted_clusters, start=1)}
    df_cluster['Cluster'] = df_cluster['Cluster'].map(cluster_mapping)
    
    # Calculer les centroïdes
    centroids_new = compute_spatial_centroids(df_cluster, cluster_col='Cluster', geometry_col='geometry')
    centroids_cluster_array = np.array([(point.x, point.y) for point in centroids_new])
    
    return df_cluster[[cle_geo, 'Cluster']]


def kmeans_with_spatial_constraint(data, spatial_coords, k,distance_mode="euclidean", lambda_penalty=0.2, n_init=30, max_iter=1000, tol=1e-4,parameter=100):
    """
    Custom k-means clustering with spatial constraints and low inertia criterion over multiple runs.

    Parameters:
    - data: numpy array, the features (e.g., PC1 to PC10).
    - spatial_coords: numpy array, the spatial coordinates (e.g., centroids or points as [x, y]).
    - k: int, number of clusters.
    - lambda_penalty: float, weight for spatial penalty.
    - n_init: int, number of times to run the k-means with different initializations.
    - max_iter: int, maximum number of iterations.
    - tol: float, convergence tolerance for centroid movement or inertia change.

    Returns:
    - clusters: array, cluster labels for each point for the best result.
    - centroids: array, final centroids in feature space for the best result.
    """
    best_inertia = np.inf
    best_labels = None
    best_centroids = None


    for _ in range(n_init):
        print(f"num {_+1}")
        # Initialize centroids randomly
        rng = np.random.default_rng()
        centroids = data[rng.choice(data.shape[0], k, replace=False)]

        # Initialize labels based on initial centroids
        labels = np.argmin(cdist(data, centroids, metric='euclidean'), axis=1)
        labels = np.random.randint(0, k, size=len(data))

        # Compute initial inertia
        inertia = np.sum(np.min(cdist(data, centroids, metric='euclidean'), axis=1))

        for iteration in range(max_iter):
            # Compute distances to centroids in feature space
            feature_distances = cdist(data, centroids, metric='euclidean')

            # Calculate the spatial_distance 
            if distance_mode == "squared":
                # Cluster variance: normalized sum of squared distances to the centroid
                spatial_centroids = np.array([
                    spatial_coords[labels == i].mean(axis=0) if (labels == i).sum() > 0 else np.zeros(2)
                    for i in range(k)
                ])
                # Calculate squared Euclidean distances between spatial coordinates and centroids
                spatial_distances = cdist(spatial_coords, spatial_centroids, metric='euclidean') ** 2
        
            elif distance_mode == "euclidean":
                spatial_centroids = np.array([
                    spatial_coords[labels == i].mean(axis=0) if (labels == i).sum() > 0 else np.zeros(2)
                    for i in range(k)
                ])
                spatial_distances = cdist(spatial_coords, spatial_centroids, metric='euclidean')
                # Return the mean of these distances as compactness

            elif distance_mode == "neighbor":
                dist_matrix = cdist(spatial_coords, spatial_coords, metric='euclidean')
            
                spatial_distances = np.zeros((len(data), k))
                for label in range(0, k):
                    # Find the indices of the rows corresponding to the current label
                    label_indices = np.where(labels == label)[0]
                        # Iterate over the rows in dist_matrix that belong to the current label
                    for idx in range(len(data)):
                        # Get the distances from the current row to all other rows
                        distances_to_others = dist_matrix[label_indices,idx]
                        # Find the indices of the 10 smallest distances (excluding itself, index `idx`)
                        closest_indices = sorted(distances_to_others)[:parameter]
                        # Store the 10 closest distances in the corresponding column of closest_distances
                        spatial_distances[idx, label - 1] = np.sum(closest_indices)
                   
            else:
                raise ValueError("Invalid compactness mode.")

            # Combine feature distances and compactness penalties
            combined_distances = feature_distances/feature_distances.mean() + lambda_penalty  * spatial_distances/spatial_distances.mean()

            # Combine feature and spatial distances
            #combined_distances = feature_distances/feature_distances.mean() + lambda_penalty  * spatial_distances/spatial_distances.mean()

            # Assign clusters based on combined distances
            new_labels = np.argmin(combined_distances, axis=1)

            # Recalculate centroids
            new_centroids = np.array([data[new_labels == i].mean(axis=0) if (new_labels == i).sum() > 0 else centroids[i]
                                      for i in range(k)])
            
            # Compute new inertia
            new_inertia = np.sum(np.min(cdist(data, new_centroids, metric='euclidean'), axis=1))

            # Check convergence (inertia improvement and centroid movement)
            if np.array_equal(new_labels, labels) or np.linalg.norm(new_centroids - centroids) < tol:
                break
            if np.abs(new_inertia - inertia) < tol:  # Inertia change is too small to improve further
                break

            labels = new_labels
            centroids = new_centroids
            inertia = new_inertia  # Update inertia

        # After each run, compare the inertia to find the best result
        if inertia < best_inertia:
            best_inertia = inertia
            best_labels = labels
            best_centroids = centroids

    return best_labels

def compute_spatial_centroids(df_cluster, cluster_col='Cluster', geometry_col='geometry'):
    """
    Compute spatial centroids for each cluster based on polygon geometry.

    Parameters:
    - df_cluster: GeoDataFrame containing polygon geometries and cluster labels.
    - cluster_col: Column name indicating cluster labels.
    - geometry_col: Column name containing polygon geometries.

    Returns:
    - centroids: GeoSeries containing the spatial centroids for each cluster.
    """
    # Ensure the DataFrame is a GeoDataFrame
    if not isinstance(df_cluster, gpd.GeoDataFrame):
        df_cluster = gpd.GeoDataFrame(df_cluster, geometry=geometry_col)
    
    # Group by cluster and compute the union of polygons within each cluster
    grouped = df_cluster.groupby(cluster_col)[geometry_col].apply(lambda g: g.unary_union)
    
    # Compute centroids of the unified polygons
    centroids = grouped.centroid
    
    return centroids

def assign_missing_clusters(carte_maille, df_cluster, centroids, cle_geo,mode="euclidean"):
    """
    Assign clusters to missing data based on spatial distances to cluster centroids.

    Parameters:
    - carte_maille: DataFrame with all cells, including 'cle_geo' and 'geometry' (spatial info).
    - df_cluster: DataFrame with PCA features and cluster labels.
    - centroids: numpy array, the final centroids of clusters in feature space.

    Returns:
    - df_cluster: DataFrame with missing cells assigned to the best cluster based on spatial proximity.
    """
    # Extract cle_geo for all cells in carte_maille
    all_cle_geo = carte_maille[cle_geo].values

    # Extract cle_geo for cells present in df_cluster
    cluster_cle_geo = df_cluster[cle_geo].values

    # Identify missing cells (those in carte_maille but not in df_cluster)
    missing_cells = [cle_geo for cle_geo in all_cle_geo if cle_geo not in cluster_cle_geo]

    # Get spatial coordinates of missing cells (assuming the geometry column is used here)
    missing_spatial_coords = carte_maille[carte_maille[cle_geo].isin(missing_cells)]['geometry'].apply(lambda x: x.centroid.coords[0]).values

    # Initialize an empty list to store the assigned cluster labels
    assigned_clusters = []

    # Calculate the spatial distances to the centroids and assign the closest centroid
    for spatial_coord in missing_spatial_coords:
        # Convert the spatial coordinates to a 2D array as cdist expects a 2D array
        spatial_coord_array = np.array([spatial_coord])

        # Compute the distances from the current missing cell to all centroids
        if mode == "euclidean":
            spatial_distances = cdist(spatial_coord_array, centroids, metric='euclidean')

        elif mode == "squared":
            spatial_distances = cdist(spatial_coord_array, centroids, metric='euclidean')**2
        else:
            raise ValueError("Invalid compactness mode.")

        
        # Assign the cluster with the minimum spatial distance
        closest_cluster_index = np.argmin(spatial_distances)+1

        # Append the assigned cluster index
        assigned_clusters.append(closest_cluster_index)

    # Check if the number of assigned clusters matches the number of missing cells
    if len(assigned_clusters) != len(missing_cells):
        raise ValueError(f"The number of assigned clusters ({len(assigned_clusters)}) does not match the number of missing cells ({len(missing_cells)})")

    # Now we need to add missing cells to df_cluster
    missing_cells_df = carte_maille[carte_maille[cle_geo].isin(missing_cells)]

    # Assign the clusters to the missing cells
    missing_cells_df['Cluster'] = assigned_clusters

    # Append the missing cells with the assigned clusters to df_cluster
    df_cluster = pd.concat([df_cluster, missing_cells_df[[cle_geo,'geometry','Cluster']]], ignore_index=True)

    
    # Add the assigned clusters to df_cluster for the missing cells
    df_cluster.loc[df_cluster[cle_geo].isin(missing_cells), 'Cluster'] = assigned_clusters

    return df_cluster

def determiner_k(df_pca,n_init=30,max_cluster=20):
    
    # Étape 2 : Déterminer le nombre optimal de clusters (facultatif)
    # Utiliser l'inertie ou l'Elbow Method (méthode du coude) pour trouver le bon k
    
    inertia = []
    k_range = range(1,max_cluster)
    for k in k_range:
        kmeans = KMeans(n_clusters=k,n_init=n_init)
        kmeans.fit(df_pca)
        inertia.append(kmeans.inertia_)
        
    plt.figure(figsize=(8, 6))
    plt.plot(k_range, inertia, marker='o')
    plt.xlabel('Nombre de clusters (k)')
    plt.ylabel('Inertie')
    plt.title('Méthode du coude pour déterminer k')
    plt.show()

def etudier_composition_cluster(df,df_cluster,dico_taxo,col_value='nombreObs',cle_ID='cdRef',cle_geo='codeMaille10Km'):

    df_cluster=df_cluster.reset_index()
    df_cluster=df_cluster[[cle_geo,'Cluster']]
    df_composition=pd.merge(df,df_cluster,on=[cle_geo],how='right')
    df_grouped=df_composition.groupby(["Cluster",cle_ID])[col_value].sum()
    df_grouped=df_grouped.reset_index()
    
    df_grouped=pd.merge(df_grouped,dico_taxo,on=cle_ID)

    # Calculer le total des observations pour chaque cluster
    df_grouped['totalObs'] = df_grouped.groupby('Cluster')[col_value].transform('sum')

    # Calculer le total des observations pour chaque cluster
    df_grouped['nombreObs_unique'] = 1
    
    # Calculer la colonne % composition
    df_grouped['% composition'] = round((df_grouped[col_value] / df_grouped['totalObs']) * 100,1)
    
    # Trier les espèces par nombre d'observations pour chaque cluster
    sorted_df = df_grouped.sort_values(['Cluster', col_value], ascending=[True, False])
    

    return sorted_df

def calculer_inertie(df_pca, labels):
    """Calcule l'inertie totale d'un clustering."""
    df_pca = df_pca.copy()  # Éviter d'écraser l'original
    labels = np.array(labels)  # Conversion en array NumPy

    unique_labels = np.unique(labels)
    inertia = 0

    for label in unique_labels:
        # Sélection sécurisée des points du cluster
        cluster_points = df_pca.loc[labels == label]  
        
        if cluster_points.empty:
            continue  # Évite les erreurs si un cluster est vide

        centroid = cluster_points.mean(axis=0)  # Calcul du centre du cluster
        inertia += ((cluster_points - centroid) ** 2).sum().sum()  # Double somme

    return inertia