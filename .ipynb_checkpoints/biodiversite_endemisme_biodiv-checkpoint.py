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
def calculer_WE(df_input,col_valeur='nombreObs',cle_geo='codeMaille10Km',cle_ID='cdRef'):
    grouped = df_input.groupby(cle_ID).agg(aire_repartition=('nombreObs_unique', 'sum'))
    grouped['WE_prep']=1/grouped['aire_repartition']
    df_input=pd.merge(df_input,grouped,on=cle_ID)
    grouped_WE = df_input.groupby(cle_geo).agg(indice_d_endemisme=('WE_prep', 'sum'))
        
    return grouped_WE

def calculer_indices(df_input,col_valeur='nombreObs',cle_geo='codeMaille10Km',cle_ID='cdRef'):
    groupes_nombre_especes=calculer_nombre_especes(df_input,cle_geo,cle_ID)
    groupes_nombre_observations=calculer_nombre_observations(df_input,cle_geo,cle_ID)
    grouped_shannon=calculer_shannon(df_input,col_valeur,cle_geo,cle_ID)
    grouped_simpson=calculer_simpson(df_input,col_valeur,cle_geo,cle_ID)
    grouped_WE=calculer_WE(df_input,col_valeur,cle_geo,cle_ID)
    df_indice=pd.merge(groupes_nombre_especes,groupes_nombre_observations,on=cle_geo)
    df_indice=pd.merge(df_indice,grouped_shannon,on=cle_geo)
    df_indice=pd.merge(df_indice,grouped_simpson,on=cle_geo)
    df_indice=pd.merge(df_indice,grouped_WE,on=cle_geo)
    return df_indice