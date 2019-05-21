from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA
import pandas as pd
import scanpy as sc
import numpy as np

def pca(adata):
    #calculate pca
    pca = PCA(n_components=50)
    result_pca = pca.fit_transform(adata.X)
    result_pca = np.array(result_pca, dtype = np.float32)
    adata.obsm['X_pca'] = result_pca
    return result_pca


def bbknn(adata, neighbors_within_batch):
    result_pca = pca(adata)
    batch_unique = list(set(adata.obs['batch']))  #['0','1']
    knn_distances = np.zeros((adata.shape[0],neighbors_within_batch*len(batch_unique)))
    knn_indices = np.copy(knn_distances).astype(int)

    for i in range(len(batch_unique)):
        #reference batch
        batch_id = batch_unique[i]   #get batch id for ref_batch
        bool_idx = adata.obs['batch'] == batch_id  #get booleen index for reference batch
        ref_batch_pca = result_pca[bool_idx]  #using booleen index to get pca data for reference batch
        ref_batch_idx = np.arange(adata.shape[0])[bool_idx] #create a booleen index for ref_batch to map back to pca matrix

        for j in range(len(batch_unique)):
            #querry batch
            batch_id = batch_unique[j] #get batch id for query_batch
            bool_idx = adata.obs['batch'] == batch_id  #get booleen index for query batch
            query_batch_pca = result_pca[bool_idx]  #using booleen index to get pca data for query batch
            query_batch_idx = np.arange(adata.shape[0])[bool_idx]  #create a booleen index for query_batch to map back to pca matrix

            D = pairwise_distances(X= query_batch_pca, Y=ref_batch_pca)  #calculate pairwise_distances between query batch and ref_batch

            #sort distances and neighbors
            neighbors = np.argsort(D, axis=1)[0:,0:neighbors_within_batch]  #get indices for n nearest neighbors
            sorted_D = np.sort(D, axis=1)[0:,0:neighbors_within_batch]  #get distance for n nearest neighbors

            for n in range(neighbors.shape[0]):
                for k in range(neighbors.shape[1]):
                    temp_neighbor = neighbors[n,k]
                    neighbors[n,k]=ref_batch_idx[temp_neighbor]  #map n nearest neighbors to pca indices

            col_range = np.arange(i*neighbors_within_batch, (i+1)*neighbors_within_batch)
            knn_indices[query_batch_idx[:,None], col_range[None,:]]= neighbors
            knn_distances[query_batch_idx[:,None], col_range[None,:]] = sorted_D

    return knn_indices, knn_distances


#batch_unique=['0','1']
neighbors_within_batch=6
knn_indices, knn_distances = bbknn(adata, neighbors_within_batch)
#calculating connectivities
knn_connectivities = sc.neighbors.compute_connectivities_umap(knn_indices, knn_distances, n_obs=15476,n_neighbors=12)
adata.uns['neighbors']['connectivities'] = knn_connectivities[1]
adata.uns['neighbors']['distances'] = knn_connectivities[0]
sc.tl.louvain(adata)
sc.tl.umap(adata)
sc.pl.umap(adata,  color=['Cell type', 'batch', 'louvain'])  #ploting
