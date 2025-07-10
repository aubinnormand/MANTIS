#affichage_carte_biodiv.py

from fonctions_annexes_biodiv import (
    round_to_sig
)
from normalisation_biodiv import normaliser_log
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import geopandas as gpd
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def afficher_fond_carte(nom, dictionnaire_cartes,source_fond='OpenStreetMap'):
    if nom not in dictionnaire_cartes:
        print(f"Erreur : '{nom}' non trouvé dans la liste des cartes.")
        return
    
    params = dictionnaire_cartes[nom]
    fig, ax1 = configurer_carte(
        source_fond,
        params["center_x"], params["center_y"], params["height"],
        zoom=params["zoom"], fig_size=(params["size_x"], params["size_y"])
    )
    
    return fig, ax1


def configurer_carte(basemap_type, center_x, center_y, height, zoom=6, fig_size=(13, 9)):
    """
    Configure une carte avec un fond de carte spécifié, centrée sur un point donné et une largeur donnée.
    
    Parameters:
        basemap_type (str): Le type de fond de carte (OpenStreetMap, CartoDB, Esri, NASAGIBS, GeoportailSatellite, etc.).
        center_x (float): Coordonnée X (longitude en projection) du centre de la carte.
        center_y (float): Coordonnée Y (latitude en projection) du centre de la carte.
        width (float): Largeur de la carte en unités de projection.
        zoom (int): Niveau de zoom pour le fond de carte.
        aspect_ratio (float): Ratio largeur/hauteur pour conserver les proportions (1.5 par défaut).
        
    Returns:
        fig, ax1: Figure et axe configurés avec le fond de carte.
    """
    # Calculer l'aspect ratio en fonction de la taille de la figure
    aspect_ratio = fig_size[0] / fig_size[1]
    
    # Calculer les limites de la carte en fonction du centre et de la largeur
    ylim_min = center_y - height / 2
    ylim_max = center_y + height / 2
    width = height * aspect_ratio  # Calculer la largeur pour conserver les proportions
    
    xlim_min = center_x - width / 2
    xlim_max = center_x + width / 2
    
    # Basemaps disponibles
    basemaps = {
        'OpenStreetMap': ctx.providers.OpenStreetMap.Mapnik,
        'OpenStreetMapTopo': "https://tile.thunderforest.com/topo/{z}/{x}/{y}.png?apikey=a73842edce474a5991e6d45b3d684470",
        'CartoDB': ctx.providers.CartoDB.Positron,
        'Esri': ctx.providers.Esri.WorldImagery,
        'NASAGIBS': ctx.providers.NASAGIBS.BlueMarble,
        'GeoportailSatellite': ctx.providers.GeoportailFrance.orthos,
        'GeoportailDepartements': ctx.providers.GeoportailFrance.parcels,
        'GeoportailRegion': ctx.providers.GeoportailFrance.Cartes_Naturalearth
    }
    

    # Créer la figure et ajuster les dimensions
    fig, ax1 = plt.subplots(figsize=fig_size)
    ax1.set_xlim(xlim_min, xlim_max)
    ax1.set_ylim(ylim_min, ylim_max)

    # Ajouter le fond de carte
    ctx.add_basemap(ax1, source=basemaps[basemap_type], zoom=zoom)
    ax1.axis('off')  # Désactiver les axes

    return fig, ax1

def ajouter_couche_SIG(fig, ax1, couche_sig,edgecolor="white",facecolor='none',
                       linewidth=0.5,linestyle='--',alpha=1,
                       with_label=False,col_label='NOM_SITE',fontsize=8,label_color="white"):
    # Charger les données géographiques des départements
    couche_sig = couche_sig.to_crs(epsg=3857)

    couche_sig.plot(ax=ax1, edgecolor=edgecolor, facecolor=facecolor, linewidth=linewidth, linestyle=linestyle,alpha=alpha)
    if with_label is True:
        couche_sig['centroid'] = couche_sig.geometry.centroid
        for idx, row in couche_sig.iterrows():
            if 'adhésion' not in row[col_label].lower():
            # Placer le label au centroïde de chaque géométrie
                ax1.text(row['centroid'].x, row['centroid'].y, row[col_label], 
                        fontsize=fontsize, ha='center', va='center', color=label_color,fontname='Arial')
    return fig, ax1

def ajouter_couche_continue(fig, ax1, df_inpn,df_geo,col_valeur,cle_geo,quantile_inf=0,quantile_sup=1, borne_min=None,borne_max=None,cmap_choice='viridis',val_alpha=0.75,missing_color=None,colorbar_choice=True,log_values=False):

    # Exemple d'ajout de données continues
    grouped = df_inpn.groupby(cle_geo)[col_valeur].sum().reset_index()

    if log_values is True:
        grouped=normaliser_log(grouped, observation_col=col_valeur,seuil_min=0)
        col_valeur=col_valeur+ '_log'
    
    # Filtrer les groupes où la somme est différente de 0
    grouped = grouped[grouped != 0].reset_index()

    # Ajouter les données géométriques
    sum_values_df=pd.merge(df_geo,grouped,on=[cle_geo],how='left').reset_index()

    gdf = gpd.GeoDataFrame(sum_values_df)
    gdf = gdf.to_crs(epsg=3857)

    # Déterminer les valeurs minimales et maximales de, l
    if borne_max is not None:
        gdf.loc[gdf[col_valeur]>borne_max,col_valeur]=borne_max
    if borne_min is not None:
        gdf.loc[gdf[col_valeur]<borne_min,col_valeur]=borne_min

    # Check for NaN values in the specified column
    if gdf[col_valeur].isnull().all():
        raise ValueError(f"The column '{col_valeur}' contains only NaN values. Cannot compute quantiles.")
    elif gdf[col_valeur].isnull().any():
        print(f"Warning: The column '{col_valeur}' contains some NaN values. These will be ignored for quantile calculations.")
    
    # Now calculate the quantiles, ignoring NaNs by default
    val_min = round_to_sig(gdf[col_valeur].quantile(quantile_inf, interpolation='linear'), direction='down')
    val_max = round_to_sig(gdf[col_valeur].quantile(quantile_sup, interpolation='linear'), direction='up')
      
    quart_val=val_min+(val_max-val_min)/4
    demi_val=val_min+(val_max-val_min)*2/4
    troisquart_val=val_min+(val_max-val_min)*3/4

    # Configurer `missing_kwds` seulement si `missing_color` est spécifié
    missing_kwds = {"color": missing_color, "label": "Valeur manquante"} if missing_color else None
    
    gdf.plot(ax=ax1,column=col_valeur, cmap=cmap_choice, legend=False, alpha=val_alpha,
             linewidth=0.1, edgecolor='gray',vmin=val_min,vmax=val_max,missing_kwds=missing_kwds)

    if colorbar_choice is True:
    # Ajouter une colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap_choice,norm=plt.Normalize(vmin=val_min, vmax=val_max))
        extend_param = 'both' if quantile_sup < 1 and quantile_inf > 0 else 'max' if quantile_sup < 1 else 'min' if quantile_inf > 0 else 'neither'
        cbar = fig.colorbar(sm, ax=ax1, shrink=0.6,ticks=[val_min,round_to_sig(quart_val, direction='nearest'),round_to_sig(demi_val, direction='nearest'),round_to_sig(troisquart_val, direction='nearest'),val_max],extend=extend_param)  # Réduire la taille de la colorbar
        cbar.ax.tick_params(labelsize=10)  # Augmenter la taille de la police de la colorbar
        cbar.set_label(col_valeur,fontsize=12)
        cbar.ax.yaxis.set_label_position('left') 
        
    return fig, ax1
    
def ajouter_couche_discrete(fig, ax1,df_cluster,df_geo,col_valeur='Cluster',cle_geo='codeMaille10Km',cmap_choice='plasma',val_alpha=0.75,missing_color=None,
                            legend_choice=True,loc_legend='upper left',pos_legend=(0.03, 0.55)):
    df=pd.merge(df_geo,df_cluster,on=[cle_geo],how='left')
    # Ajouter les données géométriques
    gdf = gpd.GeoDataFrame(df)
    gdf = gdf.to_crs(epsg=3857)
    
    # Déterminer les valeurs minimales et maximales de, l
    val_min = gdf[col_valeur].min()
    val_max = gdf[col_valeur].max()
    
    # Configurer `missing_kwds` seulement si `missing_color` est spécifié
    missing_kwds = {"color": missing_color, "label": "Valeur manquante"} if missing_color else None

    # Tracer la carte des départements en utilisant la colonne var pour la coloration
    gdf.plot(column=col_valeur, cmap=cmap_choice, legend=False, alpha=val_alpha, ax=ax1,missing_kwds=missing_kwds)

    values = df_cluster[col_valeur].unique()
    values=np.sort(values)

    values_color = [sorted(values).index(x) + 1 for x in values]
    
    # Obtenir la colormap
    cmap = plt.get_cmap(cmap_choice)
    
    # Calculer les positions des couleurs dans la colormap
    positions = np.linspace(0, 1, len(values))

    # Extraire les couleurs aux positions spécifiées
    colors = [cmap(pos) for pos in positions]

    if legend_choice is True:
        patches = []
        # Creating legend with color box 
        for value in values_color:
            patch= mpatches.Patch(color=colors[value-1], label=f'Cluster {values[value-1]}') 
            patches.append(patch)
        if missing_color is not None:
            patch_missing = mpatches.Patch(color=missing_color, label="Pas de données")
            patches.append(patch_missing)
        # Afficher la légende avec tous les patches
        ax1.legend(handles=patches,loc=loc_legend, bbox_to_anchor=pos_legend,fontsize=12)
        
    return fig, ax1

def ajouter_couche_point(fig, ax1,df_inpn,df_geo,col_valeur='nombreObs',cle_geo='codeMaille10Km',color_dot='black',size_dot=2,legend_choice=True):
  # Regrouper les données par 'codeMaille10Km' et 'geometry' et calculer la somme des observations
    grouped =df_inpn.groupby([cle_geo])[col_valeur].sum()
    grouped=grouped.reset_index()
    
    # Convertir l'objet Series résultant en DataFrame
    sum_values_df=pd.merge(df_geo,grouped,on=[cle_geo],how='right')
    sum_values_df_avec = sum_values_df.reset_index()
    gdf = gpd.GeoDataFrame(sum_values_df)
    # Choisir un système de coordonnées projetées approprié (par exemple, EPSG:3857 pour les coordonnées Web Mercator)
    gdf = gdf.to_crs(epsg=3857)

    ax1.scatter(gdf.geometry.centroid.x, gdf.geometry.centroid.y, color=color_dot, s=size_dot)

    if legend_choice is True:
        circle = Line2D([0], [0], marker='o', color='w', label='Circle',
                                markerfacecolor=color_dot, markersize=10),
        ax1.legend(handles=circle,labels=['Espèce recensée'], loc='upper right', fontsize=12)
        
    return fig, ax1

def ajouter_couche_statut(fig, ax1, df_statut, df_geo, statut_colors,col_valeur='statut',cle_geo='codeMaille10Km', val_alpha=0.75, legend_choice=True,loc_legend="upper right"):
    # Fusionner les données géométriques et le DataFrame des statuts
    df = pd.merge(df_geo, df_statut, on=[cle_geo], how='right')
    gdf = gpd.GeoDataFrame(df)
    gdf = gdf.to_crs(epsg=3857)

    # Vérifier si les statuts du DataFrame existent dans le dictionnaire
    missing_statuts = set(gdf[col_valeur].unique()) - set(statut_colors.keys())
    if missing_statuts:
        raise ValueError(f"Les statuts suivants manquent dans statut_colors: {missing_statuts}")

    # Ajouter une colonne de couleur dans le GeoDataFrame en fonction du statut
    gdf['color'] = gdf[col_valeur].map(statut_colors)

    # Tracer la carte en utilisant la couleur assignée pour chaque statut
    gdf.plot(color=gdf['color'], alpha=val_alpha, ax=ax1)

    # Créer la légende
    if legend_choice:
        patches = [
            mpatches.Patch(color=color, label=status) for status, color in statut_colors.items()
        ]
        ax1.legend(handles=patches, loc=loc_legend, fontsize=10)
        
    return fig, ax1

def afficher_carte_defaut(titre,df_filt,df_geo,couche_sig,
                          col_valeur='nombreObs',cle_geo='codeMaille10Km',quantile_inf=0,quantile_sup=1,borne_min=None,borne_max=None,
                           cmap_choice='viridis',basemap_type='GeoportailSatellite',save_path=None):
    fig, ax=configurer_carte(basemap_type)
    ajouter_couche_continue(fig, ax, df_filt,df_geo,col_valeur,cle_geo, quantile_inf,quantile_sup,borne_min,borne_max,cmap_choice)
    if couche_sig is not None:
        fig, ax=ajouter_couche_SIG(fig, ax,couche_sig)
    ax.set_title(titre, fontsize=16)  # Taille de la police définie à 16
    if save_path is not None:
        fig.savefig(save_path+'/'+titre+'.png', dpi=300, bbox_inches='tight')  # Enregistre au format PNG avec une résolution de 300 DPI
    return fig, ax