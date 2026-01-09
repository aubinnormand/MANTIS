#normalisation_biodiv.py 
# Fonctions de normalisation
import numpy as np

# Fonction pour caper les valeurs par la 10ème plus grande
def normaliser_par_maille(df, cle_geo='codeMaille10Km', observation_col='nombreObs'):
    # Calcul de la somme maximale parmi les groupes
    col_norm=observation_col + '_norm_par_maille'
    target_mean = df.groupby([cle_geo])[observation_col].sum().mean()
    
    somme_obs= df.groupby([cle_geo])[observation_col].transform('sum')

    # Appliquer le facteur de normalisation à chaque observation
    
    df[col_norm] = (df[observation_col] / somme_obs) * target_mean
    print('Il y a en moyenne '+str(target_mean.round(1))+' observations par maille')

    return df

def normaliser_par_aire(df, carte_maille, cle_geo='codeMaille10Km',
                        observation_col='nombreObs'):
    """
    Normalise les observations par aire (en km²) en utilisant directement 
    la colonne 'area_km2' de carte_maille.

    df : dataframe avec au moins [cle_geo, observation_col]
    carte_maille : GeoDataFrame avec [cle_geo, 'area_km2']
    observation_col : nom de la colonne contenant les observations
    """

    # Vérifier que l'aire est bien présente
    if 'area_km2' not in carte_maille.columns:
        raise ValueError("La colonne 'area_km2' doit être présente dans carte_maille.")

    # Fusion des données pour associer les aires aux observations
    df = df.merge(carte_maille[[cle_geo, 'area_km2']], on=cle_geo, how='left')

    # Calcul de la densité moyenne d’observations (nb obs / km²)
    densite_moyenne = (
        df.groupby(cle_geo)[observation_col].sum()
        /
        df.groupby(cle_geo)['area_km2'].first()
    ).mean()

    # Données par maille
    somme_obs = df.groupby(cle_geo)[observation_col].transform('sum')
    aire = df.groupby(cle_geo)['area_km2'].transform('first')

    # Colonne normalisée
    col_norm = observation_col + '_norm_par_aire'
    df[col_norm] = (df[observation_col] / somme_obs) * densite_moyenne * aire

    print(f"Il y a en moyenne {densite_moyenne.round(1)} observations par km².")

    return df



def normaliser_par_espece(df, code_col='cdRef', observation_col='nombreObs'):
    # Calcul de la somme maximale parmi les groupes
    col_norm=observation_col + '_norm_par_espece'
    # Calcul de la somme des nombreObs par nomScientifique
    somme_obs_par_nom = df.groupby(code_col)[observation_col].transform('sum')
    
    # Normalisation de chaque nombreObs
    df[col_norm] = (df[observation_col] / somme_obs_par_nom) * 10000

    print('Le nombre d observations total par espèce est fixé à 10 000')
    return df

def normaliser_par_clade(df, clade_col='regne', observation_col='nombreObs'):
    col_norm=observation_col + '_norm_par_'+clade_col

    target_mean = df.groupby([clade_col])[observation_col].sum().mean()
    
    somme_obs= df.groupby([clade_col])[observation_col].transform('sum')

    # Appliquer le facteur de normalisation à chaque observation
    
    df[col_norm] = (df[observation_col] / somme_obs) * target_mean
    print('Il y a en moyenne '+str(target_mean.round(1))+' observations par '+clade_col)


    return df
    
def normaliser_par_maille_et_clade(df, cle_geo='codeMaille10Km', clade_col='regne', observation_col='nombreObs'):
    col_norm=observation_col + '_norm_par_maille_et_'+clade_col
    target_mean = df.groupby([cle_geo,clade_col])[observation_col].sum().mean()
    
    somme_obs_par_maille_clade= df.groupby([cle_geo,clade_col])[observation_col].transform('sum')

    # Appliquer le facteur de normalisation à chaque observation
    
    df[col_norm] = (df[observation_col] / somme_obs_par_maille_clade) * target_mean
    print('Il y a en moyenne '+str(target_mean.round(1))+' observations par '+clade_col+' par maille')

    return df

def normaliser_par_aire_et_clade(df, carte_maille, cle_geo='codeMaille10Km',
                                 clade_col='regne', observation_col='nombreObs'):
    """
    Normalise les observations par maille et par clade en utilisant directement 
    la colonne 'area_km2' de carte_maille.

    df : dataframe avec au moins [cle_geo, clade_col, observation_col]
    carte_maille : GeoDataFrame avec [cle_geo, 'area_km2']
    clade_col : colonne représentant le clade (ex: regne, classe)
    observation_col : colonne contenant le nombre d'observations
    """

    # Vérifier que l'aire est bien présente
    if 'area_km2' not in carte_maille.columns:
        raise ValueError("La colonne 'area_km2' doit être présente dans carte_maille.")

    col_norm = f"{observation_col}_norm_par_aire_et_{clade_col}"

    # Fusion des données pour associer les aires aux observations
    df = df.merge(carte_maille[[cle_geo, 'area_km2']], on=cle_geo, how='left')

    # Somme des observations par maille et clade
    somme_obs_par_maille_clade = df.groupby([cle_geo, clade_col])[observation_col].transform('sum')
    somme_aire_par_maille = df.groupby([cle_geo])['area_km2'].transform('first')  # Même aire pour toutes les lignes d'une maille

    # Moyenne cible pondérée par l'aire
    target_mean = (
        df.groupby([cle_geo, clade_col])[observation_col].sum()
        /
        df.groupby([cle_geo])['area_km2'].first()
    ).mean()

    # Normalisation
    df[col_norm] = (df[observation_col] / somme_obs_par_maille_clade) * target_mean * somme_aire_par_maille

    # Nettoyage
    df = df.drop(columns=['area_km2'], errors='ignore')

    print(f"Il y a en moyenne {target_mean.round(1)} observations par {clade_col} et par km².")

    return df


    
def normaliser_unique(df,observation_col='nombreObs'):
    col_norm=observation_col+ '_unique'
    df[col_norm] = df[observation_col].apply(lambda x: 1 if x > 0 else 0)
    return df

def normaliser_log(df, observation_col='nombreObs',seuil_min=0):
    col_norm=observation_col+ '_log'
    min_value = df[df[observation_col] > seuil_min][observation_col].min()
    df[col_norm] = df[observation_col]/min_value
    # Appliquer la fonction de normalisation sur le DataFrame
    df[col_norm] = df[col_norm].apply(lambda x: np.log10(x) if x > seuil_min else 0)

    return df

def normaliser_par_periode(df, observation_col='nombreObs',col_norm='periode'):
    # Calcul de la somme maximale parmi les groupes
    nom_col_norm=observation_col + '_norm_par_periode'
    target_mean = df.groupby([col_norm])[observation_col].sum().mean()
    
    somme_obs= df.groupby([col_norm])[observation_col].transform('sum')

    # Appliquer le facteur de normalisation à chaque observation
    
    df[nom_col_norm] = (df[observation_col] / somme_obs) * target_mean
    print('Il y a en moyenne '+str(target_mean.round(1))+' observations par periode dans '+observation_col)

    return df
    
    