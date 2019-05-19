from sklearn.metrics import pairwise_distances

def pca(matrix):
    #calculate pca
    pca = PCA(n_components=50)
    result_pca = pca.fit_transform(matrix)
    result_pca = np.array(result_pca, dtype = np.float32)
    adata.obsm['X_pca'] = result_pca


def bbknn(batch_unique,neighbors_within_batch):
    knn_distances = np.zeros((adata.shape[0],neighbors_within_batch*len(batch_unique)))
    knn_indices = np.copy(knn_distances).astype(int)

    for i in range(len(batch_unique)):
        #reference batch
        #get batch id

        batch_id = batch_unique[i]
        bool_idx = adata.obs['batch'] == batch_id  #booleans of sixe n_obs (15467)
        ref_batch_pca = result_pca[bool_idx]
        ref_batch_idx = np.arange(adata.shape[0])[bool_idx]

        for j in range(len(batch_unique)):
            #querry batch
            batch_id = batch_unique[j]
            bool_idx = adata.obs['batch'] == batch_id  #booleans of sixe n_obs (15467)
            query_batch_pca = result_pca[bool_idx]
            query_batch_idx = np.arange(adata.shape[0])[bool_idx]

            #get distance between querry and ref batch
            #sort distances get neighbors

            D = pairwise_distances(X= query_batch_pca, Y=ref_batch_pca)

            #sort distances get neighbors
            neighbors = np.argsort(D, axis=1)[0:,0:neighbors_within_batch]
            sorted_D = np.sort(D, axis=1)[0:,0:neighbors_within_batch]

            for n in range(neighbors.shape[0]):
                for k in range(neighbors.shape[1]):
                    temp_neighbor = neighbors[n,k]
                    neighbors[n,k]=ref_batch_idx[temp_neighbor]

            col_range = np.arange(i*neighbors_within_batch, (i+1)*neighbors_within_batch)
            knn_indices[query_batch_idx[:,None], col_range[None,:]]= neighbors
            knn_distances[query_batch_idx[:,None], col_range[None,:]] = sorted_D

    return knn_indices, knn_distances

result_pca = pca(adata.X)
batch_unique=['0','1']
neighbors_within_batch=6
knn_indices, knn_distances = bbknn(batch_unique,neighbors_within_batch)
#calculating connectivities
knn_connectivities = sc.neighbors.compute_connectivities_umap(knn_indices, knn_distances, n_obs=15476,n_neighbors=12)
adata.uns['neighbors']['connectivities'] = knn_connectivities[1]
adata.uns['neighbors']['distances'] = knn_connectivities[0]
sc.tl.louvain(adata)
sc.tl.umap(adata)
sc.pl.umap(adata,  color=['Cell type', 'batch'])
