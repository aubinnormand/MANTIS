# correlation_prediction.py

# Fonctions de corrélation
import pandas as pd
import numpy as np
from fonctions_annexes_biodiv import generer_dictionnaire_taxonomie
import time
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.stats as stats

def calculer_correlation_sujet(df,var,clade,sujet,methode='pearson',cle_geo='codeMaille10Km',cle_ID='cdRef'):
    df_pivot = df.pivot_table(index=cle_geo, columns=clade, values=var, aggfunc='sum', fill_value=0)
    
    # Select the column corresponding to the species of interest
    species_data = df_pivot[sujet]
        
    # Compute correlations only between the selected species and the others
    correlation_with_species = df_pivot.corrwith(species_data,method=methode)
    
    # Sort the results to find the most correlated species
    correlated_species = correlation_with_species.sort_values(ascending=False)
    
    df_corr = pd.DataFrame(correlated_species)
    df_corr = correlated_species.reset_index()
    df_corr.columns = [clade, 'Coeff_corr']
    df_dico=generer_dictionnaire_taxonomie(df,cle_ID)
    df_corr=pd.merge(df_corr,df_dico,on=clade)

    return df_corr

def calculer_matrice_correlation(df_input,var='nombreObs',m='kendall',cle_geo='codeMaille10Km',cle_ID='cdRef'):
    pivot_df = df_input.pivot_table(index=cle_geo, columns=cle_ID, values=var, fill_value=0)
    correlation_matrix = pivot_df.corr(method=m)
  
    return correlation_matrix

def recalculer_nombreObs_par_correlation(df_input,df_corr,col_valeur,cle_sujet,with_specie=1,cle_geo='codeMaille10Km',cle_ID='cdRef'):

    # Créer une nouvelle colonne pour stocker les résultats
    colonne_resultat = f"{col_valeur}_predit"
    
    if with_specie==0:
        df_corr.loc[df_corr[cle_ID]==cle_sujet,'Coeff_corr']=0
    if with_specie==1:
        df_corr.loc[df_corr[cle_ID]==cle_sujet,'Coeff_corr']=1
    
    # Effectuer la jointure des deux dataframes sur 'nomScientifique'
    df_merged = pd.merge(df_input,df_corr[[cle_ID,'Coeff_corr']], on=cle_ID, how='outer')
    
    # Créer une nouvelle colonne en multipliant 'nombreObs' par 'coeff'
    df_merged[colonne_resultat] = df_merged[col_valeur] * df_merged['Coeff_corr']

    df_sujet_predit=df_merged.groupby(cle_geo)[colonne_resultat].sum()
    df_sujet_predit=df_sujet_predit.reset_index()

    return df_sujet_predit

def calculer_seuil(df_input,df_sujet_predit,cle_sujet,col_valeur,seuil_observation=1,n_sigma=2,cle_geo='codeMaille10Km',cle_ID='speciesKey'):
    colonne_resultat = f"{col_valeur}_predit"
    if n_sigma==0:
        seuil_prediction=0
    else:
        df_global_avec_presence = df_input[(df_input[cle_ID]==cle_sujet)&(df_input['nombreObs']>=seuil_observation)]
        liste_mailles_avec_sujet = df_global_avec_presence[cle_geo].unique()


        df_predit_avec_presence=df_sujet_predit[(df_sujet_predit[cle_geo].isin(liste_mailles_avec_sujet))]

        seuil_prediction=df_predit_avec_presence[colonne_resultat].mean()-n_sigma*df_predit_avec_presence[colonne_resultat].std()
    
    return seuil_prediction

def calculer_prediction(df_input, df_corr, colonne_valeurs, cle_geo='codeMaille10Km',cle_ID='cdRef'):
    """
    Cette fonction multiplie les valeurs d'une matrice de corrélation filtrée par une colonne spécifique du DataFrame df_input.
    Le résultat est stocké dans une nouvelle colonne nommée `colonne_valeurs_predit`.
    
    :param df_input: DataFrame d'entrée contenant les observations.
    :param df_corr: DataFrame contenant la matrice de corrélation.
    :param colonne_valeurs: Nom de la colonne dans df_input à multiplier avec la matrice de corrélation.
    :param cle_geo: Colonne indiquant les mailles dans df_input (par défaut: 'codeMaille10Km').
    :return: df_input avec une nouvelle colonne `colonne_valeurs_predit`.
    """

    df_input = df_input.sort_values(by=cle_ID,ascending=True)
    
    # Créer une nouvelle colonne pour stocker les résultats
    colonne_resultat = f"{colonne_valeurs}_predit"
    df_input[colonne_resultat] = np.nan

    # Boucler sur chaque maille unique
    total_maille = len(df_input[cle_geo].unique())  # Nombre total de mailles à traiter
    print(f"nombre total de mailles: {total_maille}")
    
    current_maille = 0  # Compteur d'avancement pour les mailles
    last_percent = -1  # Suivi du dernier pourcentage affiché

    # Boucler sur chaque maille unique
    for maille in df_input[cle_geo].unique():

        # Filtrer df_input pour la maille en cours
        df_maille = df_input[df_input[cle_geo] == maille]

        # Extraire la liste des espèces présentes dans cette maille
        liste_especes = df_maille[cle_ID].unique()

        # Filtrer la matrice de corrélation pour garder les espèces présentes dans la maille
        try:
            #df_corr_filtered = df_corr.loc[liste_especes, liste_especes]

            
            df_corr.columns = df_corr.columns.astype(int)
            df_corr.index = df_corr.columns

            common_species = df_corr.columns.intersection(liste_especes)
            df_corr_filtered = df_corr.loc[common_species, common_species]
            
            # Effectuer la multiplication matricielle
            df_input.loc[df_input[cle_geo] == maille, colonne_resultat] = np.dot(df_corr_filtered.values, df_maille[colonne_valeurs].values)
            
        except KeyError:
            # Gestion d'erreur si certaines espèces ne sont pas présentes dans la matrice de corrélation
            print(f"Espèces non trouvées dans la matrice de corrélation pour la maille {maille}")

        # Mise à jour de l'avancement

        current_maille += 1
        percent = int((current_maille / total_maille) * 100)

        # Afficher le pourcentage d'avancement uniquement lorsque le pourcentage change
        if percent > last_percent:
            print(f"Progress: {percent}%")
            last_percent = percent

    return df_input

def recherche_espece_absente(df_complete,col_values_predit,cle_geo='codeMaille10Km',cle_ID='cdRef'):
    
    df_complete_grouped=df_complete.groupby(cle_ID)[['nombreObs',col_values_predit]].sum()
    df_complete_grouped=df_complete_grouped.reset_index()
    
    df_especes_absentes=df_complete_grouped[df_complete_grouped['nombreObs']==0]
    df_especes_absentes = df_especes_absentes.sort_values(by=col_values_predit,ascending=False)
    
    df_dico=generer_dictionnaire_taxonomie(df_complete,cle_ID)
    df_especes_absentes=pd.merge(df_especes_absentes,df_dico,on=cle_ID)

def quadratic_model(x, a, b, c):
    return a * np.power(x, 2) + b * x + c

def power_model(x, a, b):
    return a * np.power(x, b)

def linear_model(x, a, b):
    return a * x + b

def exp_model(x, a, b):
    return np.exp(a * x + b)

def normal_cdf(x, a, b):
    return stats.norm.cdf(x, loc=a, scale=b)

def poisson_cdf(x, a):
    return stats.poisson.cdf(x, mu=a)

def prepare_data(df_global, df_complet_predit_espece, cle_geo, espece, col_values_corr):
    colonne_predit = f"{col_values_corr}_predit"
    df_global_avec_presence = df_global[df_global['species'] == espece]
    
    df_merge = pd.merge(df_global_avec_presence[[cle_geo, 'species', 'nombreObs', col_values_corr]],
                        df_complet_predit_espece[[cle_geo, colonne_predit]],
                        on=cle_geo, how='outer').fillna(0)
    df_merge = df_merge[df_merge['nombreObs'] >= 1]
    
    x = df_merge[col_values_corr].values
    y = df_merge[colonne_predit].values / df_merge[colonne_predit].max()
    
    return x, y

def fit_and_plot(x, y, col_values_corr, espece, loi='all'):
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, alpha=0.7, color='b', label='Données')
    x_range = np.linspace(min(x), max(x), 100)
    
    coefficients = {}
    
    if loi in ['linear', 'all']:
        params, _ = curve_fit(linear_model, x, y, maxfev=10000)
        a, b = params
        plt.plot(x_range, linear_model(x_range, a, b), color='r', label='Courbe de tendance affine')
        coefficients['linear'] = (a, b)
    
    if loi in ['power', 'all']:
        params, _ = curve_fit(power_model, x, y, maxfev=10000)
        a, b = params
        plt.plot(x_range, power_model(x_range, a, b), color='g', label='Courbe de tendance puissance')
        coefficients['power'] = (a, b)
    
    if loi in ['poisson', 'all']:
        params, _ = curve_fit(poisson_cdf, x, y, maxfev=10000)
        a = params[0]
        plt.plot(x_range, poisson_cdf(x_range, a), color='orange', label='Courbe de tendance poisson')
        coefficients['poisson'] = (a,)
    
    if loi in ['normal', 'all']:
        params, _ = curve_fit(normal_cdf, x, y, maxfev=10000)
        a, b = params
        plt.plot(x_range, normal_cdf(x_range, a, b), color='pink', label='Courbe de tendance normale')
        coefficients['normal'] = (a, b)
    
    plt.xlabel(col_values_corr)
    plt.ylabel(f"{col_values_corr}_predit")
    plt.title(f'Scatter Plot avec régressions pour {espece}')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    for key, value in coefficients.items():
        print(f'Les coefficients de la loi {key} pour {espece} sont {[round(v, 3) for v in value]}')

    return coefficients

def plot_residuals(x, y, a_lineaire, b_lineaire):
    y_pred_lineaire = linear_model(x, a_lineaire, b_lineaire)
    residuals = y - y_pred_lineaire
    
    plt.figure(figsize=(10, 6))
    plt.hist(residuals, bins=30, density=True, alpha=0.6, color='b', label="Résidus")
    
    mu, std = np.mean(residuals), np.std(residuals)
    xmin, xmax = plt.xlim()
    x_range = np.linspace(xmin, xmax, 100)
    pdf = stats.norm.pdf(x_range, mu, std)
    
    plt.plot(x_range, pdf, 'r', linewidth=2, label="Densité normale ajustée")
    plt.text(0.7 * xmax, max(pdf) * 0.8, f"μ = {round(mu, 2)}\nσ = {round(std, 2)}", 
             fontsize=12, color='red', bbox=dict(facecolor='white', alpha=0.6))
    
    plt.xlabel('Résidus')
    plt.ylabel('Densité')
    plt.title('Distribution des résidus et ajustement à une loi normale')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    print(f'Les paramètres de la loi normale ajustée aux résidus sont : mu={round(mu, 2)}, sigma={round(std, 2)}')

def inverse_power_function(y, a, b):
    return (y / a) ** (1 / b)

def inverse_linear_function(y, a, b):
    return (y - b) / a

def inverse_exp_function(y, a, b):
    return (np.log(y) - b) / a

def inverse_poisson_function(y, lambda_poisson):
    """Approximation de l'inverse de la CDF de Poisson en utilisant une normalisation"""
    return np.round(lambda_poisson + np.sqrt(lambda_poisson) * stats.norm.ppf(y)).astype(int)

def inverse_normal_function(y, mu, sigma):
    """Inverse de la CDF de la loi normale (fonction quantile)"""
    return mu + sigma * stats.norm.ppf(y)
    

def appliquer_transformation(df, colonne_predit, methode, coefficients):
    colonne_predit_modelisation = 'nombreObs_predit_glm_'+methode
    fonctions = {
        'power': inverse_power_function,
        'linear': inverse_linear_function,
        'exp': inverse_exp_function,
        'poisson': inverse_poisson_function,
        'normal': inverse_normal_function
    }
    
    if methode in fonctions:
        df[colonne_predit_modelisation] = fonctions[methode](df[colonne_predit] / df[colonne_predit].max(), *coefficients[methode])
    else:
        raise ValueError("Méthode de transformation inconnue. Choisissez parmi: power, linear, exp, poisson, normal")

    return df

    
