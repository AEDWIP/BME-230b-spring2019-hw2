from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA
import pandas as pd
import scanpy as sc
import numpy as np

def pca(adata):
    #calculate pca
    pca = PCA(n_components=50)
    adata.obsm['X_pca'] = pca.fit_transform(adata.X)


def l_k_bbknn(adata, batch_unique,neighbors_within_batch,l):

    #this method makes an l subsampling of the bbknn computed in the graph() method
    # if k=4, meaning you are finding 4 neighbors between batches, and l=2:
    ## subsample 2 neighbors from batch1 and 2 neighbors from batch2
    pca(adata) #run pca
    knn_indices, knn_distances = bbknn(adata,batch_unique,neighbors_within_batch)  #get result from bbknn 
    random_sample = random.sample(range(0,neighbors_within_batch), l)  #create indices for random sampling l from k

    if l >= neighbors_within_batch:
        raise ValueError('l cannot be equal or larger than k')
    #initiate lk_bbknn l_bbknn_indices and l_bbknn_distances matrix
    l_bbknn_indices = np.zeros((knn_indices.shape[0], len(batch_unique)*l)).astype(int)
    l_bbknn_distances = np.zeros((knn_distances.shape[0], len(batch_unique)*l))

    for i in range(len(batch_unique)):

        batch_id = batch_unique[i]  #get batch id for ref_batch
        bool_idx = adata.obs['batch'] == batch_id  #get booleen index for reference batch
        ref_batch_pca = adata.obsm['X_pca'][bool_idx]  #using booleen index to get pca data for reference batch
        ref_batch_idx = np.arange(adata.shape[0])[bool_idx] #create a booleen index for ref_batch to map back to pca matrix

        for j in range(len(batch_unique)):
            batch_id = batch_unique[j]   #get batch id for query_batch
            bool_idx = adata.obs['batch'] == batch_id  #get booleen index for query batch
            query_batch_pca = adata.obsm['X_pca'][bool_idx] #using booleen index to get pca data for query batch
            query_batch_idx = np.arange(adata.shape[0])[bool_idx]  #create a booleen index for query_batch to map back to pca matrix

            D = pairwise_distances(X=query_batch_pca, Y=ref_batch_pca)  #calculate pairwise_distances between query batch and ref_batch

            neighbors = np.argsort(D, axis=1)[0:,0:neighbors_within_batch]  #get indices for n nearest neighbors
            sorted_D = np.sort(D, axis=1)[0:,0:neighbors_within_batch]    #get distance for n nearest neighbors

            for n in range(neighbors.shape[0]):
                for k in range(neighbors.shape[1]):
                    temp_neighbor = neighbors[n,k]
                    neighbors[n,k]=ref_batch_idx[temp_neighbor]   #map n nearest neighbors to pca indices

            col_range = np.arange(i*l, (i+1)*l)
            #pass random sampled l nearest neighbors to l_bbknn_indices and distances matrix
            l_bbknn_indices[query_batch_idx[:,None], col_range[None,:]]= neighbors[:,random_sample]
            l_bbknn_distances[query_batch_idx[:,None], col_range[None,:]] = sorted_D[:,random_sample]

    return l_bbknn_indices, l_bbknn_distances

#batch_unique=['0','1']
neighbors_within_batch=6
l=3
l_bbknn_indices, l_bbknn_distances = l_k_bbknn(adata, neighbors_within_batch,l)
#calculating connectivities
lk_bbknn_connectivities = sc.neighbors.compute_connectivities_umap(l_bbknn_indices, l_bbknn_distances, n_obs=15476,n_neighbors=2*l)
adata.uns['neighbors']['connectivities'] = lk_bbknn_connectivities[1]
adata.uns['neighbors']['distances'] = lk_bbknn_connectivities[0]
sc.tl.louvain(adata)
sc.tl.umap(adata)
sc.pl.umap(adata,  color=['Cell type', 'batch','louvain'])
