# Fonctions d'exploration des données
from fonctions_annexes_biodiv import generer_dictionnaire_taxonomie
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np

def filtrer_top(df,var='nombreObs',nmin=100,cle_ID='cdRef'):
    grouped_df=df.groupby(cle_ID)[var].sum()
    grouped_df=grouped_df.reset_index()
    top_grouped_df=grouped_df[grouped_df[var]>nmin-1]
    df_filtre_top = df[df[cle_ID].isin(top_grouped_df[cle_ID])]
    print("nombre d'espèces retenues dans le df :"+str(len(df_filtre_top[cle_ID].unique()))+
          ' ('+str(round(len(df_filtre_top[cle_ID].unique())/len(df[cle_ID].unique())*100))+'%)')
    return(df_filtre_top)
    
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

def explorer_clade(df,choix,taxon,cle_ID='cdRef'):  
    liste_clade = ['all', 'all', 'regne', 'kingdom', 'classe', 'class', 'ordre', 'order', 
                   'famille', 'family', 'genre', 'genus', 'nomScientifique', 'species', 
                   'nomVernaculaire']
    
    # Vérifier si 'choix' est présent dans la liste des clades
    if choix not in liste_clade:
        # Si 'choix' n'est pas trouvé, définir 'clade_inf' à 'regne'/'kingdom'
        print(f"Clade '{choix}' not found. Defaulting to 'regne'/'kingdom'.")
        clade_inf = 'kingdom'  # Changer pour 'regne' ou 'kingdom' comme clade par défaut
        df_filt = df
    else:
        # Trouver l'index de l'élément choisi
        index_choix = liste_clade.index(choix)
        # Trouver l'élément suivant (clade inférieur)
        if index_choix < len(liste_clade) - 1:
            clade_inf = liste_clade[index_choix + 2]  # On saute deux indices pour accéder à l'inférieur
            # Filtrer le DataFrame selon le clade et le taxon choisi
            df_filt = df[df[choix] == taxon]
        else:
            clade_inf = None  # Pas d'élément suivant si 'choix' est le dernier élément dans la liste
        
    # Extraire les taxons inférieurs uniques et les observations
    if clade_inf:
        taxons_inf = df_filt[clade_inf].unique()
    else:
        taxons_inf = []  # Aucun taxon inférieur si clade_inf est None
    
    # Calculer le nombre d'observations et le nombre d'espèces uniques
    n_obs = int(df_filt['nombreObs'].sum())
    n_especes = len(df_filt[cle_ID].unique())
    # Afficher avec des séparateurs de milliers
    n_especes_formatte = "{:,}".format(n_especes)
    #print("Nombre d'espèces:", n_especes_formatte)
    
    # Afficher avec des séparateurs de milliers
    n_Obs_formatte = "{:,}".format(n_obs)
    #print("Nombre d'observations:", n_Obs_formatte)
    
    #Creer un dataframe avec comme colonnes : le nom des taxons inférieurs, le nombre d'espèces dans ce taxon et le nombre d'observation
    df_explo = pd.DataFrame()
    df_explo['taxons']=taxons_inf
    
    nobs = df_filt.groupby([clade_inf])['nombreObs'].sum()
    
    # Grouper par famille et compter les noms scientifiques uniques
    nesp = df_filt.groupby(clade_inf)[cle_ID].nunique().reset_index()
    nesp.columns = [clade_inf, 'nombreEspèces']
    
    # Renommer la colonne pour plus de clarté
    resultat=pd.merge(nesp,nobs,on=clade_inf)
    resultat['Ratio Obs/Esp']=round(resultat['nombreObs']/resultat['nombreEspèces'],1)
    # Trier le DataFrame par la colonne nombreObs
    resultat = resultat.sort_values(by='nombreEspèces', ascending=False).reset_index(drop=True)
    
    # Taille de la figure
    fig, ax1 = plt.subplots(figsize=(14, 8))
    
    # Positions des barres sur l'axe x
    indices = np.arange(len(resultat[clade_inf]))
    largeur_barres = 0.4  # Largeur des barres
    
    # Création du premier axe pour les noms scientifiques uniques
    ax1.set_yscale('log')
    bar1 = ax1.bar(indices - largeur_barres/2, resultat['nombreEspèces'], largeur_barres, label='nombreEspèces', color='b')
    ax1.set_xlabel(clade_inf)
    ax1.set_ylabel('nombreEspèces', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    
    # Créer un deuxième axe pour la somme des observations, partageant le même axe x
    ax2 = ax1.twinx()
    ax2.set_yscale('log')
    bar2 = ax2.bar(indices + largeur_barres/2, resultat['nombreObs'], largeur_barres, label='nombreObs', color='r')
    ax2.set_ylabel('nombreObs', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    
    # Ajuster les ticks de l'axe x pour qu'ils correspondent aux familles
    ax1.set_xticks(indices)
    ax1.set_xticklabels(resultat[clade_inf], rotation=45, ha='right')
    
    # Ajouter une légende combinée
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    
    # Ajuster la disposition pour éviter le chevauchement des labels
    #plt.title(f"Nombre de noms scientifiques uniques et somme des observations par famille pour l'ordre {taxon}")
    plt.tight_layout()
    plt.savefig('figure.png')
    # Afficher le graphique
    plt.show()

    # Afficher le DataFrame trié
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
    