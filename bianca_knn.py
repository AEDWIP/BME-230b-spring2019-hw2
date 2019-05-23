from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA
import pandas as pd
import scanpy as sc
import numpy as np


def pca(adata):
    #calculate pca
    pca = PCA(n_components=50)
    adata.obsm['X_pca'] = pca.fit_transform(adata.X)


def get_distances(matrix):
    D = pairwise_distances(X=matrix, Y=matrix)
    return D

def get_neighbors(D,k):
    knn_indices = np.argsort(D, axis=1)[0:,1:k+1]
    knn_distances = np.sort(D, axis=1)[0:,1:k+1]
    return knn_indices, knn_distances

def get_umap_connectivities(knn_i, knn_d, k):
    knn_connectivities = sc.neighbors.compute_connectivities_umap(knn_i, knn_d, n_obs=15476,n_neighbors=k)
    connectivities = knn_connectivities[1]
    distances = knn_connectivities[0]
    return connectivities, distances

def knn(adata, k):
    adata.uns['neighbors']['connectivities'] = None
    adata.uns['neighbors']['distances'] = None
    pca(adata)
    D = get_distances(adata.obsm['X_pca'])
    knn_indices, knn_distances = get_neighbors(D,k)
    connectivities, distances = get_umap_connectivities(knn_indices, knn_distances, k)
    #pass connectivities, distances to scanpy
    adata.uns['neighbors']['connectivities'] = connectivities
    adata.uns['neighbors']['distances'] = distances
    sc.tl.umap(adata)
    sc.pl.umap(adata,  color=['Cell type', 'batch'])

adata = sc.read('/Users/biancaxue/Documents/class/BME230B/HW2/PBMC.merged.h5ad',
                delimiter='\t',cache=True)
#adata.var_names_make_unique()
k=15
knn(adata, k)
