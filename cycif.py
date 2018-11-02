# -*- coding: utf-8 -*-
# Author Yunguan 'Jake' Wang 
# 07/16/2018
"""
Cy-cif segmented psudeo single-cell data Cell State Calling script
"""
import os
import pandas as pd
import numpy as np
import yaml
# import matplotlib
# # Force matplotlib to not use any Xwindows backend.
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.cluster import KMeans
from hdbscan import HDBSCAN
import umap
import argparse

def avg_fc_gen(dfdata, df_ctype, 
               celllist=None, 
               class_col=0):
    """
    Calculate class average and ratio of class averages
    """
    data = dfdata.copy()
    classtype_sub = df_ctype.copy()
    if celllist is None:
        data = data.loc[classtype_sub.index]
    else:
        data = data.loc[celllist]
        classtype_sub = classtype_sub.loc[celllist]
    data['group'] = classtype_sub.iloc[:, class_col].values
    cavg = data.groupby('group').mean().transpose()
    tar_cell_pops = classtype_sub.iloc[:, class_col].unique()
    df_foldchange = pd.DataFrame(index = cavg.index,columns=cavg.columns)
    for pop in tar_cell_pops:
        if cavg.shape[1] == 1:
            df_foldchange = cavg
        else:
            cells = classtype_sub[classtype_sub.iloc[:, class_col] != pop].index.tolist()
            overall_avg = data.loc[cells].mean()  
            df_foldchange[pop] = cavg[pop].values - overall_avg.values
    return cavg, df_foldchange

def get_data(filepath):
    """
    read data into a dataframe, the input can be a local directory, a dropbox link, or a synapse id.
    """
    filepath = './input/' + filepath
    if filepath[-3:] == 'csv':
        data = pd.read_csv(filepath)
    elif filepath[-3:] == 'txt':
        data = pd.read_table(filepath)
    elif filepath[-4:] == 'xlsx':
        data = pd.read_excel(filepath)
    else:
        raise FilenameError('File name not supported!')
    return data


def preprocess_and_normalize(df_raw_expr, keep_location_data=False):
    new_index = ['cell_'+str(x) for x in df_raw_expr.index]
    df_raw_expr.index = new_index
    # valid_cols = [x for x in df_raw_expr.columns if ('hoechst' not in x.lower())&(x.lower() not in ['x','y','frame', 'area','circ','xt','yt','rows','cols'])]
    # downsampling the data by location but taking 
    if keep_location_data:
        df_location = df_raw_expr[['Xt','Yt']]
    else:
        df_location = None
    # df_raw_expr = df_raw_expr[valid_cols]
    df_raw_expr = np.log10(df_raw_expr+1)
    # Use quantile transformation to mitigate outlier effects
    # ?df_norm = preprocessing.quantile_transform(df_raw_expr,n_quantiles=1000,copy=True,output_distribution='normal')
    df_norm = preprocessing.robust_scale(df_raw_expr,quantile_range=(1, 99),copy=True)
    df_norm = pd.DataFrame(df_norm,index=df_raw_expr.index,columns=df_raw_expr.columns)
    return df_raw_expr, df_location, df_norm

def dimention_reduction(df_norm):
    """
    Dimention reduction using the Uniform Manifold Approximation and Projection method.
    n_neighbors: This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default.
    min_dist: This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are more evenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5, with 0.1 being a reasonable default.
    metric: This determines the choice of metric used to measure distance in the input space. A wide variety of metrics are already coded, and a user defined function can be passed as long as it has been JITd by numba.
    
    By default using 5 components, 15 neighbors and min_dist at 0.01. 
    https://github.com/lmcinnes/umap;  
    """
    # To do: make these parameters tunable. 
    umap_embedding = umap.UMAP(n_components=5,n_neighbors=15, min_dist=0.01, metric='correlation')
    df_umap = umap_embedding.fit_transform(df_norm)
    return df_umap

def clustering(df_low_dim,
                algorithm='KMeans', 
                K_means_n_clusters = 15,
                hdbscan_min_cluster_size = 20):
    """
    Perform clustering on dimention reduced data. 
    Clustering algorithms can be KMeans or HDBSCAN.
    
    """
    if algorithm == 'HDBSCAN':
        clustering = HDBSCAN(min_cluster_size=hdbscan_min_cluster_size)
    else:
        clustering = KMeans(K_means_n_clusters)
    labels = clustering.fit_predict(df_low_dim)
    df_labels = pd.Series(['Cluster ' + str(x) for x in labels])
    return df_labels

def plot_clustering(df_low_dim, df_labels):
    """
    Make scatter plots of dimention reduced data with cluster designation
    """
    plt.ioff()
    for cluster in sorted(df_labels.unique()):
        cells_idx = df_labels[df_labels==cluster].index.values
        plt.scatter(df_low_dim[cells_idx,0],df_low_dim[cells_idx,1],s = 0.01,label=cluster)
    plt.legend(markerscale=50,bbox_to_anchor=(1, 0.9))
    plt.savefig('Clustering_on_2D.png',bbox_inches='tight')

def plot_expr_on_2D(df_2d,df_raw_expr):
    """
    Make a grid of scatter plots of original data projected to a 2D space determined by UMAP.
    Each marker is visualized in a subplot and each point in the subplot is colored based on its original expression. 
    """
    plt.ioff()
    nrows = int(np.ceil(len(df_raw_expr.columns)/5))
    fig, axes = plt.subplots(nrows, 5, sharex='col', sharey='row',squeeze=False,figsize=(25,5*nrows))
    axes = axes.ravel()
    axes_idx = 0
    for col in df_raw_expr.columns:
        subplot = axes[axes_idx].scatter(df_2d[:,0],df_2d[:,1],c = df_raw_expr[col].values,s = 0.01, cmap='bwr',label = col)
        axes[axes_idx].legend(frameon=False)
        fig.colorbar(subplot,ax=[axes[axes_idx],axes[axes_idx]],pad=-0.05,extend='both',format='%.1f')
    # for i in range(nrows):
    #     for j in range(5):
    #         colname = df_raw_expr.columns[5*i+j]
    #         subplot = ax[i,j].scatter(df_2d[:,0],df_2d[:,1],c = df_raw_expr[colname].values,s = 0.01, cmap='bwr',label = colname)
    #         ax[i,j].legend(frameon=False)
    #         fig.colorbar(subplot,ax=[ax[i,j],ax[i,j]],pad=-0.05,extend='both',format='%.1f')
    plt.savefig('raw expr on 2D.png',bbox_inches='tight')

def get_top_markers(df_raw_expr,df_labels):
    """
    Differential expression analysis based on fold change against population only. 
    Very naive currently, may be improved upon request.
    Get top 5 markers and theri fold change.
    """
    ctype = pd.DataFrame(df_labels.values,index = df_raw_expr.index, columns = ['cluster'])
    ctype.to_csv('Cluster assignments.csv')
    _, df_fc = avg_fc_gen(df_raw_expr,ctype)
    df_report = pd.DataFrame(index = df_fc.columns,columns = ['Top 5 markers', 'Markers exp'])
    for x in df_report.index:
        top_5_markers = df_fc[x].sort_values(ascending = False).index[:5]
        top_5_markers_exp = df_fc.transpose().loc[x,top_5_markers]
        df_report.loc[x,'Top 5 markers'] = ', '.join(top_5_markers)
        df_report.loc[x,'Markers exp'] = ', '.join([str(np.round(x,2)) for x in top_5_markers_exp])
    df_report.to_csv('Top markers per cluster.csv')

if __name__ == '__main__':
    # Parse all necessary arguments
    if not os.path.exists('output'):
        os.makedirs('output')
    with open('./config/cycif.yml') as f:
        config = yaml.load(f)

    file_path = config['file_path']
    algorithm = config['algo']
    n_clusters = config['nc']
    mix_cluster_size = config['mcs']
    cfg_genes = config['markers']

    parser = argparse.ArgumentParser(description='Cy-cif segmented psudeo single-cell data Cell State Calling script. ')
    parser.add_argument('file_path', nargs='?',default=file_path, help='Target data matrix of samples by features. Take segmented cycif data and output clustering figures, expression figures and top markers. Will use example file in /input folder by default.')
    parser.add_argument('-algo', metavar='Clustering algorithm', nargs='?', default=algorithm, help="Clustering algorithm, select from ['KMeans','HDBSCAN']. KMeans by default.")
    parser.add_argument('-mcs', metavar='Minimal cluster size', nargs='?', default=mix_cluster_size, type=int, help='Minimal cluster size for HDBSCAN clustering. 15 by default.')
    parser.add_argument('-nc', metavar='Number of clusters', nargs='?', default=n_clusters, type=int, help='Number of clusters for KMeans clustering. 15 by default.')
    args = parser.parse_args()
    # Read data from dropbox, synpase or local
    df_raw = get_data(args.file_path)
    cm_genes = [x for x in cfg_genes if x in df_raw.columns]
    if len(cm_genes) < len(cfg_genes):
        print('Dropped genes in config file only: {}'.format(', '.join(genes_in_cfg_only)))
        genes_in_cfg_only = [x for x in cfg_genes if x not in df_raw.columns]
    elif len(cm_genes) < df_raw.shape[1]:
        genes_in_input_only =[x for x in df_raw.columns if x not in cfg_genes]
        print('Dropped genes in input file only: {}'.format(', '.join(genes_in_input_only)))

    df_raw =df_raw[cm_genes]
    # preprocessing
    df_raw, _, df_norm = preprocess_and_normalize(df_raw)
    # Dimention reduction.
    df_umap = dimention_reduction(df_norm)
    # Clustering based on selected algorithm.
    df_labels = clustering(df_umap,
                            algorithm=args.algo, 
                            K_means_n_clusters=args.nc,
                            hdbscan_min_cluster_size=args.mcs)
    # Generating plots and reports.
    os.chdir('output')
    plot_clustering(df_umap,df_labels)
    plot_expr_on_2D(df_umap,df_raw)
    get_top_markers(df_raw,df_labels)
    