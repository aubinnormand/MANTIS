# Fonctions d'exploration des données
from fonctions_annexes_biodiv import generer_dictionnaire_taxonomie
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np

def filtrer_top_global(df, var='nombreObs', nmin=100, cle_ID='cdRef', display=True):
    grouped_df = df.groupby(cle_ID)[var].sum().reset_index()
    top_grouped_df = grouped_df[grouped_df[var] >= nmin]
    df_filtre_top = df[df[cle_ID].isin(top_grouped_df[cle_ID])]
    
    if display:
        nb = len(df_filtre_top[cle_ID].unique())
        total = len(df[cle_ID].unique())
        print(f"[global] {nb} espèces retenues ({round(nb/total*100)}% du total)")
    
    return df_filtre_top


def filtrer_top_mailles(df, var='nombreObs', ntop=100, cle_ID='cdRef', cle_geo='codeMaille10km', display=True):
    # Classement par maille et espèces
    top_species = (
        df.groupby([cle_geo, cle_ID])[var].sum()
        .reset_index()
        .sort_values([cle_geo, var], ascending=[True, False])
    )
    # Récupérer les ntop par maille
    top_species = top_species.groupby(cle_geo).head(ntop)
    
    # Liste unique des espèces
    especes_conservees = top_species[cle_ID].unique().tolist()
    
    # Filtrer le DataFrame initial
    df_filt = df[df[cle_ID].isin(especes_conservees)].copy()
    
    if display:
        nb = len(df_filt[cle_ID].unique())
        total = len(df[cle_ID].unique())
        print(f"[mailles] {nb} espèces retenues ({round(nb/total*100)}% du total)")
    
    return df_filt


def filtrer_top(df, method='global',var='nombreObs',var_global=None,var_mailles=None, nmin_global=100, ntop_mailles=10, cle_ID='cdRef', cle_geo='codeMaille10km', display=True):
    if var_global is None:
        var_global=var
    if var_mailles is None:
        var_mailles=var
        
    if method == 'global':
        df_filt = filtrer_top_global(df, var=var_global, nmin=nmin_global, cle_ID=cle_ID, display=display)

    elif method == 'mailles':
        df_filt = filtrer_top_mailles(df, var=var_mailles, ntop=ntop_mailles, cle_ID=cle_ID, cle_geo=cle_geo, display=display)

    elif method == 'combined':
        # filtrage global
        df_global = filtrer_top_global(df, var=var_global, nmin=nmin_global, cle_ID=cle_ID, display=False)
        # filtrage par mailles
        df_mailles = filtrer_top_mailles(df, var=var_mailles, ntop=ntop_mailles, cle_ID=cle_ID, cle_geo=cle_geo, display=False)
        
        # Récupérer les ensembles d'espèces
        set_mailles = set(df_mailles[cle_ID].unique())
        set_global = set(df_global[cle_ID].unique())
        
        # Intersection = espèces communes
        liste_espece = set_mailles | set_global

        df_filt = df[df[cle_ID].isin(liste_espece)].copy()
        
        if display:
            nb = len(liste_espece)
            total = len(df[cle_ID].unique())
            print(f"[combined] {nb} espèces retenues ({round(nb/total*100)}% du total)")

    else:
        print("Choisir méthode global, mailles ou combined")
        return None

    return df_filt


def afficher_top_especes(df,df_dico,col_values='nombreObs',cle_ID='cdRef'):
    # Sum each taxon column individually
    grouped = df.groupby([cle_ID])[col_values].sum()
    grouped=grouped.reset_index()
    grouped=pd.merge(grouped,df_dico,on=cle_ID)
    # Sort the individual taxon sums in descending orde
    grouped = grouped.sort_values(by=col_values,ascending=False)
    grouped = grouped.reset_index(drop=True)
    # Afficher le DataFrame trié
    return grouped
    
def chercher_espece(df, dico_taxo, var, col_valeur='nombreObs', cle_ID='cdRef'):
    # Appliquer un masque pour détecter la présence de 'var' dans n'importe quelle colonne de type chaîne de caractères
    mask = df.map(lambda x: isinstance(x, str) and var.lower() in x.lower() if isinstance(x, str) else False).any(axis=1)
    
    # Filtrer les lignes correspondant à la recherche
    recherche_df = df[mask]
    
    # Grouper par cle_ID et sommer les valeurs de col_valeur
    grouped = recherche_df.groupby(cle_ID)[col_valeur].sum()
    grouped_df = grouped.reset_index()
    
    # Fusionner avec le dictionnaire taxonomique
    grouped_df = pd.merge(grouped_df, dico_taxo, on=cle_ID)
    
    # Trier et réinitialiser l'index
    grouped_df = grouped_df.sort_values(by=col_valeur, ascending=False).reset_index(drop=True)
    
    return grouped_df

def explorer_clade(df, choix, taxon, cle_ID='cdRef', colorder='nombreEspèces'):  
    liste_clade = ['all', 'all', 'regne', 'kingdom', 'embranchement','phylum','classe', 'class', 'ordre', 'order', 
                   'famille', 'family', 'genre', 'genus', 'nomScientifique', 'species', 
                   'nomVernaculaire']
    
    # Vérifier si 'choix' est présent dans la liste des clades
    if choix not in liste_clade:
        print(f"Clade '{choix}' not found. Defaulting to 'regne'/'kingdom'.")
        clade_inf = 'kingdom'
        df_filt = df
    else:
        index_choix = liste_clade.index(choix)
        if index_choix < len(liste_clade) - 1:
            clade_inf = liste_clade[index_choix + 2]  # On saute deux indices
            df_filt = df[df[choix] == taxon]
        else:
            clade_inf = None
        
    # Extraire les taxons inférieurs uniques
    if clade_inf:
        taxons_inf = df_filt[clade_inf].unique()
    else:
        taxons_inf = []
    
    # Calcul des totaux
    n_obs = int(df_filt['nombreObs'].sum())
    n_especes = len(df_filt[cle_ID].unique())
    
    print(f"Nombre d'espèces : {n_especes:,}")
    print(f"Nombre d'observations : {n_obs:,}")
    
    # DataFrame avec les taxons inférieurs
    df_explo = pd.DataFrame({'taxons': taxons_inf})
    
    # Inclure les NaN avec dropna=False
    nobs = df_filt.groupby([clade_inf], dropna=False)['nombreObs'].sum()
    nesp = df_filt.groupby([clade_inf], dropna=False)[cle_ID].nunique().reset_index()
    nesp.columns = [clade_inf, 'nombreEspèces']
    
    resultat = pd.merge(nesp, nobs, on=clade_inf, how='outer')
    resultat['Ratio Obs/Esp'] = round(resultat['nombreObs'] / resultat['nombreEspèces'], 1)
    
    # Remplacer les NaN du clade inférieur par un libellé explicite
    resultat[clade_inf] = resultat[clade_inf].fillna('Non renseigné')
    
    # Trier le DataFrame
    resultat = resultat.sort_values(by=colorder, ascending=False).reset_index(drop=True)
    
    # --------- Graphique ---------
    fig, ax1 = plt.subplots(figsize=(14, 8))
    
    indices = np.arange(len(resultat[clade_inf]))
    largeur_barres = 0.4
    
    # Axe 1 : nombre d'espèces
    ax1.set_yscale('log')
    bar1 = ax1.bar(indices - largeur_barres/2, resultat['nombreEspèces'], largeur_barres, label='nombreEspèces', color='b')
    ax1.set_xlabel(clade_inf)
    ax1.set_ylabel('nombreEspèces', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    
    # Axe 2 : nombre d'observations
    ax2 = ax1.twinx()
    ax2.set_yscale('log')
    bar2 = ax2.bar(indices + largeur_barres/2, resultat['nombreObs'], largeur_barres, label='nombreObs', color='r')
    ax2.set_ylabel('nombreObs', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    
    # Ticks X
    ax1.set_xticks(indices)
    ax1.set_xticklabels(resultat[clade_inf], rotation=45, ha='right')
    
    # Légendes combinées
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc='upper right')
    
    plt.tight_layout()
    plt.savefig('figure.png')
    plt.show()
    
    return resultat


def chercher_especes_protegees(df_inpn,liste_codes,cle_geo='codeMaille10Km',cle_ID='cdRef',col_valeur='nombreObs'):
    df_local=df_inpn[(df_inpn[cle_geo].isin(liste_codes))]
    df_local_especes_protegee=df_local[df_local['especeProtegee']=='True']
    grouped_local_especes_protegee=df_local_especes_protegee.groupby([cle_ID])[col_valeur].sum()
    grouped_local_especes_protegee=grouped_local_especes_protegee.reset_index()
    grouped_local_especes_protegee = grouped_local_especes_protegee.sort_values(by=col_valeur,ascending=False)
    df_dico=generer_dictionnaire_taxonomie(df_local,cle_ID)
    grouped_local_especes_protegee=pd.merge(grouped_local_especes_protegee,df_dico,on=cle_ID)
    return grouped_local_especes_protegee
    