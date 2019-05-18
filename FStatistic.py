from sklearn.decomposition import PCA
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.spatial import distance_matrix
from aedknn import knnG

k = 12
adata = sc.read('/Users/jcasaletto/PycharmProjects/BME230B/HW2/BME-230b-spring2019-hw2/PBMC.merged.h5ad')

def PCA_init(matrix):
    pca = PCA(n_components=50)
    result_pca = pca.fit_transform(matrix)
    df = pd.DataFrame(result_pca)
    return df


def bbknn(k):  # k is the number of neighbors within each batch
    pca_total = PCA_init(adata.X)

    adata_0 = pca_total[0:8098]  # batch=0, 8098 cells
    adata_1 = pca_total[8098:]  # batch=1, 7378 cells
    batch_unique = [0, 1]

    total_data = [adata_0, adata_1]

    total_neighbors_d = []
    total_neighbors_i = []
    for i in range(len(batch_unique)):
        for j in range(len(batch_unique)):
            neighbors_d = []
            neighbors_i = []
            dist_matrix = pd.DataFrame(distance_matrix(total_data[j], total_data[i]), index=total_data[j].index,
                                       columns=total_data[i].index)
            for cell_index in range(len(dist_matrix.index)):
                temp_d = list(dist_matrix.iloc[cell_index].sort_values())[1:k + 1]
                temp_i = list(np.argsort(dist_matrix.iloc[cell_index]))[1:k + 1]
                neighbors_d.append(temp_d)
                neighbors_i.append(temp_i)
            neighbors_d = np.array(neighbors_d)
            neighbors_i = np.array(neighbors_i)
            total_neighbors_d.append(neighbors_d)  # this is the list with 4 parts
            total_neighbors_i.append(neighbors_i)

    knn_distances = np.zeros((adata.shape[0], k * len(batch_unique)))
    knn_indices = np.copy(knn_distances).astype(int)

    knn_distances[0:8098, 0:k] = total_neighbors_d[0]
    knn_distances[8098:, 0:k] = total_neighbors_d[1]
    knn_distances[0:8098, k:] = total_neighbors_d[2]
    knn_distances[8098:, k:] = total_neighbors_d[3]

    knn_indices[0:8098, 0:k] = total_neighbors_i[0]
    knn_indices[8098:, 0:k] = total_neighbors_i[1]
    knn_indices[0:8098, k:] = total_neighbors_i[2]
    knn_indices[8098:, k:] = total_neighbors_i[3]

    return knn_distances, knn_indices

knn_distances, knn_indices = bbknn(k)

knn_connectivities = sc.neighbors.compute_connectivities_umap(knn_indices, knn_distances, n_obs=15476,n_neighbors=k)

adata.uns['neighbors']['connectivities'] = knn_connectivities[1]
adata.uns['neighbors']['distances'] = knn_connectivities[0]
sc.tl.louvain(adata)  #update louvain
sc.tl.umap(adata)     #update umap
#sc.pl.umap(adata, color=['louvain'])
#sc.pl.umap(adata, color=['Cell type'])
#sc.pl.umap(adata, color=['Method'])

'''
# create list of umap points for each louvain cluster
louvainClusterPointMap = {k: [] for k in range(len(adata.obs['louvain'].cat.categories))}
for i in range(len(adata.obs)):
    louvainClusterPointMap[int(adata.obs.iloc[i]['louvain'])].append(adata.obsm['X_umap'][i])

# calculate centroids of louvain clusters
louvainClusterCentroids = []
for i in range(len(adata.obs['louvain'].cat.categories)):
    louvainClusterCentroids.append(tuple(map(lambda y: sum(y) / float(len(y)), zip(*louvainClusterPointMap[i]))))
'''


def calculate_MSE_C(adata):
    # create map of cell types to points
    cellTypePointMap = {k: [] for k in adata.obs['Cell type'].cat.categories}
    for i in range(len(adata.obs)):
        cellTypePointMap[adata.obs.iloc[i]['anno']].append(adata.obsm['X_umap'][i])

    # calculate centroids of cell types
    cellTypeCentroids = []
    for i in adata.obs['Cell type'].cat.categories:
        cellTypeCentroids.append(tuple(map(lambda y: sum(y) / float(len(y)), zip(*cellTypePointMap[i]))))

    # create map of cell types
    cellTypes = adata.obs['Cell type'].cat.categories
    cellTypeMap = {k:cellTypes[k] for k in range(len(cellTypes))}

    # calculate SSE across each cell type
    SSE_C = 0.0
    for cell_type_1 in range(len(cellTypes) - 1):
        for cell_type_2 in range(cell_type_1 + 1, len(cellTypes)):
            for coord in range(2):
                SSE_C +=(cellTypeCentroids[cell_type_1][coord] - cellTypeCentroids[cell_type_2][coord]) ** 2

    T = len(cellTypes)
    return SSE_C/(T*(T+1))


def calculate_MSE_B(adata):
    # create map of batches to points
    batchPointMap = {k: [] for k in adata.obs['batch'].cat.categories}
    for i in range(len(adata.obs)):
        batchPointMap[adata.obs.iloc[i]['batch']].append(adata.obsm['X_umap'][i])


    # create centroids of batches
    batchCentroids = {batch: tuple() for batch in adata.obs['batch'].cat.categories}
    for batch in adata.obs['batch'].cat.categories:
        batchCentroids[batch] = tuple(map(lambda y: sum(y) / float(len(y)), zip(*batchPointMap[batch])))

    SSE_B = {batch: 0.0 for batch in adata.obs['batch'].cat.categories}
    for batch in adata.obs['batch'].cat.categories:
        for point in batchPointMap[batch]:
            for coord in range(2):
                SSE_B[batch] += (point[coord] - batchCentroids[batch][coord]) ** 2

    MSE_B = 0.0
    for batch in adata.obs['batch'].cat.categories:
        MSE_B += SSE_B[batch] / len(batchPointMap[batch])

    return 1 / 4 * MSE_B

MSE_B_bb = calculate_MSE_B(adata)
MSE_C_bb = calculate_MSE_C(adata)
# now finally calculate F statistic as ratio MSE_C / MSE_B
F_stat_bb = MSE_C_bb / MSE_B_bb
print(str(F_stat_bb))


pbmc = sc.read('/Users/jcasaletto/PycharmProjects/BME230B/HW2/BME-230b-spring2019-hw2/PBMC.merged.h5ad')
myKNNG = knnG(pbmc, n_neighbors=12)
sc.tl.louvain(myKNNG.adata, resolution=1, flavor='igraph', directed=False, use_weights=True)
sc.pl.umap(myKNNG.adata, color=["louvain"])
MSE_B = calculate_MSE_B(myKNNG.adata)
MSE_C = calculate_MSE_C(myKNNG.adata)
F_stat = MSE_C / MSE_B
print(str(F_stat))


