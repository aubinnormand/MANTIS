# fonctions_annexes_biodiv.py
import pandas as pd
import math
import numpy as np
import geopandas as gpd
import re
from IPython.display import display, HTML
import matplotlib.pyplot as plt

def round_to_sig(num, sig=3, direction='nevarest'):
    if num is None:
        raise ValueError("La valeur à arrondir ne doit pas être None.")
    if num == 0:
        return 0
    # Déterminer l'ordre de grandeur du nombre
    order_of_magnitude = math.floor(math.log10(abs(num)))
    # Calculer le facteur de multiplication
    factor = 10**(sig - order_of_magnitude - 1)
    if direction == 'nearest':
        return round(num * factor) / factor
    elif direction == 'up':
        return math.ceil(num * factor) / factor
    elif direction == 'down':
        return math.floor(num * factor) / factor
    else:
        raise ValueError("Direction must be 'nearest', 'up', or 'down'")

def generer_dictionnaire_taxonomie(df,cle_ID='cdRef'):
    # Définir les colonnes attendues
    colonnes_attendues = ['speciesKey','cdRef','cdNom','taxonID','speciesID',
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

def afficher_dataframe(df, liste_colonnes="all", col_sort=None, ascending=False,n_rows=None,n_dec=1):
    # Identifier les colonnes contenant "nombreObs"
    colonnes_obs = [col for col in df.columns if 'nombreObs' in col.lower()]
    
    # Arrondir ces colonnes à 1 décimale si elles ne sont pas déjà des entiers
    for col in colonnes_obs:
        if col in df.columns and not np.issubdtype(df[col].dtype, np.integer):
            #df[col] = df[col].apply(lambda x: round(x, n_dec) if pd.notnull(x) else x)
            df[col] = df[col].round(n_dec)
            
    # Sélectionner les colonnes demandées (ou toutes si "all")
    if liste_colonnes != "all":
        df = df[liste_colonnes]

    # Supprimer les colonnes dupliquées
    df = df.loc[:, ~df.columns.duplicated()]

    # Trier si une colonne de tri est spécifiée
    if col_sort is not None:
        df = df.sort_values([col_sort], ascending=ascending)

    # Réinitialiser l'index
    df = df.reset_index(drop=True)

    # Limiter le nombre de lignes si spécifié
    if n_rows is not None:
        df = df.head(n_rows)

    return df
    
def completer_df(df_local,df_global,cle_geo='codeMaille10Km',cle_ID='cdRef'):
    # Étape 1: Récupérer les valeurs uniques de codeMaille10Km et nomScientifique
    unique_codeMaille = df_local[cle_geo].unique()
    unique_nomScientifique = df_global[cle_ID].unique()
    
    # Étape 2: Créer toutes les combinaisons possibles
    combinations = pd.MultiIndex.from_product([unique_codeMaille, unique_nomScientifique], names=[cle_geo, cle_ID]).to_frame(index=False)

    liste_nombres=[col for col in df_local.columns if 'nombreObs' in col]
    
    df_local_reduit=df_local[[cle_geo,cle_ID]+liste_nombres]
    
    # Étape 3: Fusionner avec le DataFrame original pour identifier les lignes manquantes
    df_local_complet = combinations.merge(df_local_reduit, on=[cle_geo, cle_ID], how='left')

    df_local_complet=df_local_complet.fillna(0)
    
    # Étape 4: Remplir les valeurs manquantes dans les colonnes normalisées
    df_dico=generer_dictionnaire_taxonomie(df_global,cle_ID=cle_ID)
    df_local_complet=pd.merge(df_local_complet,df_dico,on=cle_ID)
    
    return df_local_complet

def lister_mailles_dans_site(df_geo,parc_gpd,nom_parc,taux_min=0.5,cle_geo='codeMaille10Km',cle_nom_site='NOM_SITE',methode='exact'):
    # Assurons-nous que les deux GeoDataFrames sont dans le même système de coordonnées
    df_geo = df_geo.to_crs(parc_gpd.crs)
    if methode == 'contains':
        parc_gpd_filt=parc_gpd[parc_gpd[cle_nom_site].str.contains(nom_parc,case=False)]
    elif methode == 'exact':
        parc_gpd_filt=parc_gpd[parc_gpd[cle_nom_site]==nom_parc]
    else: 
        print("La méthode n'est pas valide, methode= 'contains' ou 'exact'")
        parc_gpd_filt=parc_gpd[parc_gpd[cle_nom_site]==nom_parc]

    #print(f"La zone choisiz est {parc_gpd_filt[cle_nom_site].iloc[0]}")
    
    # Calculer l'union de toutes les géométries de PNR_gpd avec union_all() au lieu de unary_union
    pnr_union = parc_gpd_filt.geometry.union_all()
    joined = gpd.sjoin(df_geo, parc_gpd, how="inner", predicate="intersects")
    # Calculer l'intersection entre chaque géométrie de carte_maille et PNR_gpd
    joined["intersection"] = joined.apply(lambda row: row["geometry"].intersection(pnr_union), axis=1)
    
    # Calculer la surface des géométries et de leur intersection
    joined["area_original"] = joined["geometry"].area
    joined["area_intersection"] = joined["intersection"].area
    
    # Filtrer pour garder les géométries dont l'intersection couvre au moins 50% de leur surface
    filtered = joined[joined["area_intersection"] >= taux_min* joined["area_original"]]
    
    # Supprimer les colonnes temporaires si elles ne sont pas nécessaires
    filtered = filtered.drop(columns=["intersection", "area_original", "area_intersection", "index_right"])
    liste_mailles=filtered[cle_geo].reset_index(drop=True)
    return liste_mailles

def ajouter_nom_site_df(df_inpn,df_geo,parc_gpd,taux_min=0.5,cle_geo='codeMaille10Km',cle_nom_site='NOM_SITE',display=True,methode='all'):
    liste_noms=parc_gpd[cle_nom_site].str.replace(r"\s*\[aire d'adhésion\]", "", regex=True, flags=re.IGNORECASE).drop_duplicates().reset_index(drop=True)
    for nom_site in liste_noms:
        liste_mailles=lister_mailles_dans_site(df_geo,parc_gpd,nom_site,taux_min,cle_geo,cle_nom_site,methode)
        if display is True:
            print(f"{len(liste_mailles)} mailles ont au moins {round(taux_min*100)}% de leur surface comprise dans le PN {nom_site}")
        df_inpn.loc[df_inpn[cle_geo].isin(liste_mailles),cle_nom_site]=nom_site
    return df_inpn


def filtrer_grille(grid, lat_min, lat_max, lon_min, lon_max):
    # Exclure les cellules hors des bornes données
    grid = grid[(
        (grid.geometry.centroid.y >= lat_min) &
        (grid.geometry.centroid.y <= lat_max) &
        (grid.geometry.centroid.x >= lon_min) &
        (grid.geometry.centroid.x <= lon_max)
    )]
    return grid

def afficher_carte_monde(world_terrestre,ecart=15):
    # Créer la figure et l'axe
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Tracer les données géographiques
    world_terrestre.plot(ax=ax, edgecolor='black', facecolor='none', alpha=0.7)
    
    # Ajouter des méridiens et des parallèles à intervalles réguliers
    # Définir les intervalles
    meridians = np.arange(-180, 181, ecart)  # Tous les 30 degrés
    parallels = np.arange(-90, 91, ecart)   # Tous les 30 degrés
    
    # Tracer les méridiens
    for meridian in meridians:
        ax.plot([meridian, meridian], [-90, 90], color='gray', linestyle='--', linewidth=0.5)
    
    # Tracer les parallèles
    for parallel in parallels:
        ax.plot([-180, 180], [parallel, parallel], color='gray', linestyle='--', linewidth=0.5)
    
    # Ajouter des étiquettes pour les méridiens et les parallèles
    for meridian in meridians:
        ax.text(meridian, -95, f"{meridian}°", color='black', ha='center', fontsize=8)
    for parallel in parallels:
        ax.text(-185, parallel, f"{parallel}°", color='black', va='center', fontsize=8)
    
    # Ajouter un titre
    ax.set_title('Carte du monde avec méridiens et parallèles', fontsize=14)
    
    # Ajuster les limites
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    
    # Afficher la carte
    plt.show()

def filtrer_geo_coord(df_biodiv, carte_maille, cle_geo,
                           lon_min=-11, lon_max=25, lat_min=32, lat_max=55):
    """
    Filtre et agrège les données de biodiversité en fonction des limites géographiques et des mailles valides.

    Paramètres :
    ------------
    - df_biodiv : DataFrame contenant les observations de biodiversité
    - carte_maille : DataFrame contenant les mailles géographiques
    - dico_taxo : DataFrame contenant les informations taxonomiques
    - cle_geo : str, nom de la colonne contenant l'identifiant des mailles
    - cle_ID : str, nom de la colonne contenant l'identifiant des espèces
    - filtrage_grille : bool, True pour appliquer un filtre géographique sur la grille
    - lon_min, lon_max, lat_min, lat_max : float, limites géographiques du filtrage

    Retourne :
    ----------
    - df_biodiv_periode : DataFrame filtré et agrégé par période
    - df_biodiv_sansperiode : DataFrame filtré et agrégé sans période
    """
    
    # Vérifier si le GeoDataFrame a un CRS projeté, sinon le reprojeter
    if carte_maille.crs is None or carte_maille.crs.is_geographic:
        print("🔄 Reprojection des données en EPSG:3857 pour un bon filtrage")
        carte_maille = carte_maille.to_crs(epsg=3857)

    # Calcul des centroides après reprojection
    carte_maille["centroid"] = carte_maille.geometry.centroid
    carte_maille_filt = carte_maille[
        (carte_maille["centroid"].y >= lat_min) & (carte_maille["centroid"].y <= lat_max) &
        (carte_maille["centroid"].x >= lon_min) & (carte_maille["centroid"].x <= lon_max)
    ]
               
    print("Filtrage géographique activé")
    carte_maille_filt = filtrer_grille(carte_maille, lat_min, lat_max, lon_min, lon_max)
    df_biodiv = df_biodiv[df_biodiv[cle_geo].isin(carte_maille_filt[cle_geo])]

    return df_biodiv

def filtrer_categorie(df, filtre, filtres):
    """
    Applique un filtre spécifique sur un DataFrame et retourne un DataFrame filtré.

    :param df: DataFrame sur lequel appliquer le filtre.
    :param filtre: Nom du filtre à appliquer (doit être une clé du dictionnaire `filtres`).
    :param filtres: Dictionnaire des filtres (modifiable si besoin).
    :return: DataFrame filtré.
    """
    # Vérifier si le filtre demandé existe
    if filtre not in filtres:
        raise ValueError(f"Filtre '{filtre}' non reconnu. Filtres disponibles : {list(filtres.keys())}")

    # Construire les conditions de filtrage
    conditions = filtres[filtre]
    filtres_list = [(df[col].isin(valeurs)) for col, valeurs in conditions.items()]

    # Appliquer le filtre (garde les lignes qui correspondent à AU MOINS une condition)
    df_filt = df[pd.concat(filtres_list, axis=1).any(axis=1)].copy()

    return df_filt

def filtrer_geographie_zone(df_raw_import, carte_maille,
                       type_site=None, site=None, 
                       cle_geo="codeMaille10Km", cle_nom_site="nom", 
                       taux_min=0.01):
    """
    Retourne la liste des mailles correspondant à une sélection géographique.
    """
    # Reprojection en Lambert 93
    type_site = type_site.to_crs(epsg=2154)
    carte_maille = carte_maille.to_crs(epsg=2154)


    liste_mailles = []

    if site:  # Vérifie qu'il y a bien une liste de sites
        for s in site:
            mailles = lister_mailles_dans_site(
                carte_maille, type_site, s,
                taux_min=taux_min, cle_geo=cle_geo, cle_nom_site=cle_nom_site, methode='exact'
            )
            liste_mailles.extend(mailles)

        # Suppression des doublons
        liste_mailles = list(set(liste_mailles))

    return liste_mailles
    



