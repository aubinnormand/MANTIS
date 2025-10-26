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

def normaliser_par_aire(df, carte_maille, cle_geo='codeMaille10Km',size_grid=None, observation_col='nombreObs'):
    # Calcul de l'aire des mailles
    carte_maille = carte_maille.to_crs(epsg=4326)
    carte_maille['aire'] = carte_maille.geometry.area
    
    # Si size_grid est défini, normaliser l'aire pour qu'elle soit à la taille de la grille
    if size_grid is not None:
        max_aire = carte_maille['aire'].max()  # Trouver la valeur maximale de l'aire
        carte_maille['aire'] = (carte_maille['aire'] / max_aire) * (size_grid ** 2)
        print(f"Les aires sont normalisées à {size_grid} km2.")
        
    # Fusion des données pour associer les aires aux observations
    df = df.merge(carte_maille[[cle_geo, 'aire']], on=cle_geo, how='left')
    
    # Calcul de la somme des observations pondérée par l'aire
    col_norm = observation_col + '_norm_par_maille'
    target_mean = (df.groupby([cle_geo])[observation_col].sum() / df.groupby([cle_geo])['aire'].first()).mean()
    
    somme_obs = df.groupby([cle_geo])[observation_col].transform('sum')
    somme_aire = df.groupby([cle_geo])['aire'].transform('first')  # Même aire pour toutes les lignes d'une maille
    
    if size_grid is not None:
        df[col_norm] = (df[observation_col] / somme_obs) * target_mean * somme_aire / (size_grid ** 2)
        df = df.drop(columns=['aire'])
        print(f"Il y a en moyenne {target_mean.round(1)} observations par km2.")
    else:
        df[col_norm] = (df[observation_col] / somme_obs) * target_mean * somme_aire
        df = df.drop(columns=['aire'])
        print(f"Il y a en moyenne {target_mean.round(1)} observations par unité d'aire.")
    
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

def normaliser_par_aire_et_clade(df, carte_maille, cle_geo='codeMaille10Km', size_grid=None,clade_col='regne', observation_col='nombreObs'):
    col_norm = f"{observation_col}_norm_par_maille_et_{clade_col}"
    
    # Calcul de l'aire des mailles
    carte_maille = carte_maille.to_crs(epsg=4326)
    carte_maille['aire'] = carte_maille.geometry.area

    # Si size_grid est défini, normaliser l'aire pour qu'elle soit à la taille de la grille
    if size_grid is not None:
        max_aire = carte_maille['aire'].max()  # Trouver la valeur maximale de l'aire
        carte_maille['aire'] = (carte_maille['aire'] / max_aire) * (size_grid ** 2)
        print(f"Les aires sont normalisées à {size_grid} km2.")
        
    # Fusion des données pour associer les aires aux observations
    df = df.merge(carte_maille[[cle_geo, 'aire']], on=cle_geo, how='left')
    
    # Calcul du nombre d'observations normalisé par maille et clade
    somme_obs_par_maille_clade = df.groupby([cle_geo, clade_col])[observation_col].transform('sum')
    somme_aire_par_maille = df.groupby([cle_geo])['aire'].transform('first')  # Même aire pour une maille donnée
    
    # Moyenne cible pondérée par l'aire
    target_mean = (df.groupby([cle_geo, clade_col])[observation_col].sum() / df.groupby([cle_geo])['aire'].first()).mean()
    
    # Normalisation en tenant compte de l'aire
    if size_grid is not None:
        df[col_norm] = (df[observation_col] / somme_obs_par_maille_clade) * target_mean * somme_aire_par_maille / (size_grid ** 2)
        df = df.drop(columns=['aire'])
        print(f"Il y a en moyenne {target_mean.round(1)} observations par {clade_col} et par km2.")
    else:
        df[col_norm] = (df[observation_col] / somme_obs_par_maille_clade) * target_mean * somme_aire_par_maille
        df = df.drop(columns=['aire'])
        print(f"Il y a en moyenne {target_mean.round(1)} observations par {clade_col} et par unité d'aire.")

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
    
    