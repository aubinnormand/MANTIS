# biodiversite_endemisme_biodiv.py

# Fonctions Indices de Biodiversité et Endémisme

import numpy as np
import pandas as pd

def calculer_nombre_especes(df_input,cle_geo='codeMaille10Km',cle_ID='cdRef'):
    df_input=df_input.drop_duplicates(subset=[cle_ID, cle_geo], keep='first')
    grouped = df_input.groupby(cle_geo)['nombreObs_unique'].agg(nombre_especes='sum').reset_index()
    return grouped

def calculer_nombre_observations(df_input,cle_geo='codeMaille10Km',cle_ID='cdRef'):
    grouped = df_input.groupby(cle_geo)['nombreObs'].agg(nombre_observations='sum').reset_index()
    return grouped

def calculer_shannon(df_input,col_valeur='nombreObs',cle_geo='codeMaille10Km',cle_ID='cdRef'):
    #Mesure de la biodiversité
    df_input = df_input.groupby([cle_geo,cle_ID])[col_valeur].sum().reset_index()
    grouped = df_input.groupby(cle_geo).agg(somme_obs=(col_valeur, 'sum'))
    df_input=pd.merge(df_input,grouped,on=cle_geo)
    df_input['prop_Obs']=df_input[col_valeur]/df_input['somme_obs']
    
    #Indice de Shannon
    df_input['Shannon_prep']=df_input['prop_Obs']*np.log10(df_input['prop_Obs'])
    grouped_shannon = df_input.groupby(cle_geo).agg(indice_de_Shannon=('Shannon_prep', 'sum'))
    grouped_shannon['indice_de_Shannon']=-grouped_shannon['indice_de_Shannon']

    return grouped_shannon

def calculer_simpson(df_input,col_valeur='nombreObs',cle_geo='codeMaille10Km',cle_ID='cdRef'):
    #Mesure de la biodiversité
    df_input = df_input.groupby([cle_geo,cle_ID])[col_valeur].sum().reset_index()
    grouped = df_input.groupby(cle_geo).agg(somme_obs=(col_valeur, 'sum'))
    df_input=pd.merge(df_input,grouped,on=cle_geo)
    df_input['prop_Obs']=df_input[col_valeur]/df_input['somme_obs']
    
    #Indice de Simpson
    df_input['Simpson_prep']=df_input['prop_Obs']**2
    grouped_simpson = df_input.groupby(cle_geo).agg(indice_de_Simpson=('Simpson_prep', 'sum'))
    grouped_simpson['indice_de_Simpson']=1-grouped_simpson['indice_de_Simpson']
    
    return grouped_simpson
    
#Mesure de l'endémicité
#Indice de Weighted Endemism (WE)
def calculer_WE(df_input, carte_maille, col_valeur='nombreObs', cle_geo='codeMaille10Km', cle_ID='speciesID'):
    """
    Calcule l'indice d'endémisme pondéré (WE) en tenant compte de l'aire réelle des cellules.
    
    Parameters
    ----------
    df_input : pd.DataFrame
        DataFrame avec au moins les colonnes [cle_geo, cle_ID].
    carte_maille : pd.DataFrame
        DataFrame avec les colonnes [cle_geo, 'area_km2'].
    col_valeur : str
        Nom de la colonne correspondant au nombre d'observations (non utilisée directement ici mais conservée).
    cle_geo : str
        Nom de la colonne identifiant les cellules.
    cle_ID : str
        Nom de la colonne identifiant les espèces.
        
    Returns
    -------
    pd.DataFrame
        DataFrame avec l'indice d'endémisme par cellule.
    """
    # Ajouter l'aire des cellules
    df_input = pd.merge(df_input, carte_maille[[cle_geo, 'area_km2']], on=cle_geo, how='left')
    
    # Calculer l'aire totale de répartition par espèce
    grouped = df_input.groupby(cle_ID).agg(aire_repartition=('area_km2', 'sum'))
    
    # Calculer WE pour chaque espèce
    grouped['WE_prep'] = 1 / grouped['aire_repartition']
    
    # Fusionner l'indice WE avec le DataFrame
    df_input = pd.merge(df_input, grouped[['WE_prep']], on=cle_ID, how='left')
    
    # Somme des WE par cellule
    grouped_WE = df_input.groupby(cle_geo).agg(indice_d_endemisme=('WE_prep', 'sum'))
    
    return grouped_WE


def calculer_margalev(df_input,col_valeur='nombreObs',cle_geo='codeMaille10Km',cle_ID='cdRef'):
    # On agrège pour avoir le nombre total d'observations par maille
    N = df_input.groupby(cle_geo)[col_valeur].sum().reset_index(name="N")
    
    # On calcule le nombre d'espèces uniques (richesse) par maille
    S = df_input.groupby(cle_geo)[cle_ID].nunique().reset_index(name="S")
    
    # On fusionne pour avoir N et S dans le même dataframe
    df_Mg = pd.merge(N, S, on=cle_geo)

    df_Mg["indice_de_Margalef"] = (df_Mg["S"] - 1) / np.log(df_Mg["N"])
    return df_Mg[[cle_geo, "indice_de_Margalef"]]   # 👈 toujours renvoyer cle_geo + R_Mg

def calculer_equitabilite_simpson(df_input, grouped_simpson, cle_geo='codeMaille10Km', cle_ID='cdRef'):
    # Richesse spécifique (nombre d'espèces uniques par maille)
    S = df_input.groupby(cle_geo)[cle_ID].nunique().reset_index(name="S")

    # Fusion avec Simpson
    df_eqS = pd.merge(grouped_simpson, S, on=cle_geo)

    # Calcul équitabilité de Simpson (E1/D)
    df_eqS["indice_equitabilite_simpson"] = (1 - df_eqS["indice_de_Simpson"]) / (1 - 1/df_eqS["S"])
    df_eqS["indice_equitabilite_simpson"]=1-df_eqS["indice_equitabilite_simpson"]

    return df_eqS[[cle_geo, "indice_equitabilite_simpson"]]

def calculer_equitabilite_heip(df_input, grouped_shannon, cle_geo='codeMaille10Km', cle_ID='cdRef'):
    # Richesse spécifique (nombre d'espèces uniques par maille)
    S = df_input.groupby(cle_geo)[cle_ID].nunique().reset_index(name="S")

    # Fusion avec Shannon
    df_eqH = pd.merge(grouped_shannon, S, on=cle_geo)

    # Calcul équitabilité de Heip
    df_eqH["indice_equitabilite_heip"] = (np.exp(df_eqH["indice_de_Shannon"]) - 1) / (df_eqH["S"] - 1)
    df_eqH["indice_equitabilite_heip"]=1-df_eqH["indice_equitabilite_heip"]

    return df_eqH[[cle_geo, "indice_equitabilite_heip"]]

def calculer_berger_parker(df_input, col_valeur='nombreObs', cle_geo='codeMaille10Km', cle_ID='cdRef', n_top=1):
    """
    Indice de Berger-Parker pour les n_top espèces les plus abondantes par maille.
    
    Parameters:
        df_input : DataFrame
        col_valeur : colonne contenant les abondances
        cle_geo : colonne de regroupement (maille)
        cle_ID : colonne identifiant les espèces
        n_top : nombre d'espèces les plus abondantes à considérer
    """
    # Nombre total d’individus par maille
    N = df_input.groupby(cle_geo)[col_valeur].sum().reset_index(name="N")
    
    # Filtrer les n_top espèces les plus abondantes par maille
    df_top = df_input.groupby(cle_geo).apply(lambda x: x.nlargest(n_top, col_valeur)).reset_index(drop=True)
    
    # Nombre total d’individus par maille (uniquement sur les n_top espèces)
    N_top = df_top.groupby(cle_geo)[col_valeur].sum().reset_index(name="N_top")
    
    # Fusion
    df_BP = pd.merge(N_top, N , on=cle_geo)
    
    # Calcul indice de Berger-Parker
    df_BP["indice_de_BergerParker"] = 1-(df_BP["N_top"] / df_BP["N"])
    
    return df_BP[[cle_geo, "indice_de_BergerParker"]]


def calculer_entropie_quadratique_taxo(df_input, col_valeur='nombreObs', 
                                                  cle_geo='codeMaille10Km', cle_ID='speciesKey'):
    """
    Calcul vectorisé de l'entropie quadratique taxonomique (Rao Q) pour toutes les mailles.
    La distance entre espèces est basée sur la taxonomie :
        - 0 si même espèce
        - 1 si même genre
        - 2 si même famille
        - 3 si même ordre
        - 4 si même classe
        - 5 si même kingdom
        - 6 sinon
    """
    df_filt = df_input.copy()
    
    # Proportion par espèce dans la maille
    df_filt['N_total'] = df_filt.groupby(cle_geo)[col_valeur].transform('sum')
    df_filt['p'] = df_filt[col_valeur] / df_filt['N_total']
    
    # Produit cartésien des espèces par maille (i,j)
    df_pairs = pd.merge(
        df_filt,
        df_filt,
        on=cle_geo,
        suffixes=('_i', '_j')
    )
    
    # Calcul vectorisé des distances taxonomiques
    def distance_taxo(row):
        if row[f'{cle_ID}_i'] == row[f'{cle_ID}_j']:
            return 0
        elif row['genus_i'] == row['genus_j']:
            return 1
        elif row['family_i'] == row['family_j']:
            return 2
        elif row['order_i'] == row['order_j']:
            return 3
        elif row['class_i'] == row['class_j']:
            return 4
        elif row['kingdom_i'] == row['kingdom_j']:
            return 5
        else:
            return 6

    df_pairs['d'] = df_pairs.apply(distance_taxo, axis=1)
    
    # Produit des proportions
    df_pairs['pq'] = df_pairs['p_i'] * df_pairs['p_j']
    df_pairs['pq_d'] = df_pairs['pq'] * df_pairs['d']
    
    # Somme par maille → entropie quadratique taxonomique
    df_rao = df_pairs.groupby(cle_geo)['pq_d'].sum().reset_index(name='entropie_quadratique_taxo')
    
    return df_rao




def calculer_indices(df_input,carte_maille, col_valeur='nombreObs', cle_geo='codeMaille10Km', cle_ID='cdRef'):
    groupes_nombre_especes = calculer_nombre_especes(df_input, cle_geo, cle_ID)
    groupes_nombre_observations = calculer_nombre_observations(df_input, cle_geo, cle_ID)
    grouped_shannon = calculer_shannon(df_input, col_valeur, cle_geo, cle_ID)
    grouped_simpson = calculer_simpson(df_input, col_valeur, cle_geo, cle_ID)
    grouped_WE = calculer_WE(df_input,carte_maille, col_valeur, cle_geo, cle_ID)
    grouped_margalev = calculer_margalev(df_input, col_valeur, cle_geo, cle_ID)
    grouped_eqSimpson = calculer_equitabilite_simpson(df_input, grouped_simpson, cle_geo, cle_ID)
    grouped_eqHeip = calculer_equitabilite_heip(df_input, grouped_shannon, cle_geo, cle_ID)
    grouped_BP = calculer_berger_parker(df_input, col_valeur, cle_geo, cle_ID)

    df_indice = pd.merge(groupes_nombre_especes, groupes_nombre_observations, on=cle_geo)
    df_indice = pd.merge(df_indice, grouped_shannon, on=cle_geo)
    df_indice = pd.merge(df_indice, grouped_simpson, on=cle_geo)
    df_indice = pd.merge(df_indice, grouped_WE, on=cle_geo)
    df_indice = pd.merge(df_indice, grouped_margalev, on=cle_geo)
    df_indice = pd.merge(df_indice, grouped_eqSimpson, on=cle_geo)
    df_indice = pd.merge(df_indice, grouped_eqHeip, on=cle_geo)
    df_indice = pd.merge(df_indice, grouped_BP, on=cle_geo)


    return df_indice


def calcul_indices_par_zone(df_zone, df_indices, cle_geo, cols, critere="moy"):
    """
    Calcule moyennes, sommes, quantiles et rangs par zone pour les colonnes spécifiées.

    Parameters
    ----------
    df_zone : pd.DataFrame
        DataFrame des zones avec colonnes 'nom' et 'liste_mailles'.
    df_indices : pd.DataFrame
        DataFrame contenant les indices par maille.
    cle_geo : str
        Nom de la colonne identifiant la maille dans df_indices.
    cols : list
        Liste des colonnes à utiliser pour calculs (moyenne, somme, etc.).
    critere : str
        Critère pour le calcul des rangs : "_moy", "_q5" ou "_q9".

    Returns
    -------
    df_indices_par_zone : pd.DataFrame
        DataFrame contenant moyennes, sommes, quantiles et rangs par zone.
    liste_colonnes : list
        Liste des colonnes à afficher (utile pour affichage ultérieur).
    """
    
    # Colonnes pour lesquelles on calcule la somme et la moyenne
    cols_somme = ["nombre_observations"]
    cols_moyenne = [c for c in cols if c not in cols_somme]
    critere=f"_{critere}"

    resultats = []

    for _, row in df_zone.iterrows():
        nom_zone = row['nom']
        mailles = row["liste_mailles"]

        # Filtrer df_indices sur les mailles de cette zone
        subset = df_indices[df_indices[cle_geo].isin(mailles)]

        if not subset.empty:
            moyennes = subset[cols_moyenne].mean()
            sommes = subset[cols_somme].sum()
            
            # Quantiles
            quantiles = subset[cols_moyenne].quantile([0.5, 0.9])  # médiane et 9ème décile
            q5 = quantiles.loc[0.5]
            q9 = quantiles.loc[0.9]

            # Stocker les résultats
            d = {"nom": nom_zone}
            d.update({f"{c}_moy": moyennes[c] for c in cols_moyenne})
            d.update({f"{c}_sum": sommes[c] for c in cols_somme})
            d.update({f"{c}_q5": q5[c] for c in cols_moyenne})
            d.update({f"{c}_q9": q9[c] for c in cols_moyenne})
            resultats.append(d)

    # Transformer en DataFrame
    df_indices_par_zone = pd.DataFrame(resultats)

    # Colonnes à classer
    cols_rang = [c for c in cols_moyenne if "indice" in c or "entropie" in c]

    # Calcul des rangs selon le critere choisi
    for c in cols_rang:
        col_a_ranger = f"{c}{critere}"
        df_indices_par_zone[f"rang_{c}"] = df_indices_par_zone[col_a_ranger].rank(
            method="min", ascending=False
        ).astype(int)

    # Rang moyen
    df_indices_par_zone["rang_moyen"] = df_indices_par_zone[[f"rang_{c}" for c in cols_rang]].mean(axis=1).astype(int)

    # Tri par rang moyen
    df_indices_par_zone = df_indices_par_zone.sort_values("rang_moyen")

    # Colonnes pour affichage
    liste_colonnes = (
        ["nom"] +
        [f"{c}_sum" for c in cols_somme] +
        [f"{c}{critere}" for c in cols_moyenne] +
        [f"rang_{c}" for c in cols_rang] +
        ["rang_moyen"]
    )

    return df_indices_par_zone, liste_colonnes

    import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

def comparer_indices(df_pred, df_obs, cle_geo, cols_indices, 
                     valeur_na=0, max_val=1e10, 
                     norm=True, method="sub"):
    """
    Compare les indices prédit et observé avec option de normalisation 
    et choix de la méthode de comparaison.
    
    Parameters
    ----------
    df_pred : pd.DataFrame
        DataFrame des valeurs prédites.
    df_obs : pd.DataFrame
        DataFrame des valeurs observées.
    cle_geo : str
        Nom de la colonne clé géographique.
    cols_indices : list of str
        Colonnes d'indices à comparer.
    valeur_na : float, default=0
        Valeur de remplacement pour les NaN après nettoyage.
    max_val : float, default=1e10
        Seuil pour filtrer les valeurs trop grandes.
    norm : bool, default=True
        Si True, applique une normalisation MinMax sur chaque DataFrame.
    method : {"sub", "div"}, default="sub"
        Méthode de comparaison :
        - "sub" : soustraction (predit - obs)
        - "div" : rapport (predit / obs)
    
    Returns
    -------
    df_compare : pd.DataFrame
        DataFrame avec colonnes comparées selon la méthode choisie.
    """
    
    def clean_df(df, cols):
        df_clean = df.copy()
        for col in cols:
            df_clean[col] = df_clean[col].replace([np.inf, -np.inf], np.nan)
            df_clean[col] = df_clean[col].apply(lambda x: x if abs(x) < max_val else np.nan)
            df_clean[col] = df_clean[col].fillna(valeur_na)
        return df_clean

    # Nettoyer
    df_pred_clean = clean_df(df_pred, cols_indices)
    df_obs_clean = clean_df(df_obs, cols_indices)
    
    # Normaliser si demandé
    if norm:
        scaler_pred = MinMaxScaler()
        df_pred_clean[cols_indices] = scaler_pred.fit_transform(df_pred_clean[cols_indices])
        
        scaler_obs = MinMaxScaler()
        df_obs_clean[cols_indices] = scaler_obs.fit_transform(df_obs_clean[cols_indices])
    
    # Initialiser df_compare
    df_compare = pd.DataFrame()
    df_compare[cle_geo] = df_pred_clean[cle_geo].values
    
    # Aligner df_obs_clean sur df_pred_clean
    df_obs_clean = df_obs_clean.set_index(cle_geo).reindex(df_compare[cle_geo]).reset_index()
    
    # Ajouter colonnes comparées
    for col in cols_indices:
        col_pred = df_pred_clean[col].values
        col_obs = df_obs_clean[col].values
        
        if method == "sub":
            df_compare[f'{col}'] = col_pred - col_obs
        elif method == "div":
            # éviter les divisions par 0
            df_compare[f'{col}'] = np.where(col_obs != 0, col_pred / col_obs, np.nan)
        else:
            raise ValueError("method doit être 'sub' ou 'div'")
    
    return df_compare


