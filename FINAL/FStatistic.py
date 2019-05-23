import logging
import scanpy as sc
from euclid_knn import KnnG
from euclid_bbknn import bbknn_graph
from bblknn import bblknn_graph
import sys


class FStatistic():
    '''
    Output: F-statistic to evaluate clustering

    methods:
        __init__(self,adata)

        calculate_MSE_C(): calculate mean-squared error across cell types

        calculate_MSE_B(): calculate mean-squared error across batches
    '''

    logger = logging.getLogger(__name__)

    ######################################################################
    def __init__(self, adata):
        '''
        input:
            adata (AnnData object)
        '''
        self._adata = adata

    def calculate_MSE_C(self):

        # create map of cell types to points
        cellTypePointMap = {k: [] for k in self._adata.obs['Cell type'].cat.categories}
        for i in range(len(self._adata.obs)):
            cellTypePointMap[self._adata.obs.iloc[i]['anno']].append(self._adata.obsm['X_umap'][i])

        # create list of cell types
        cellTypes = self._adata.obs['Cell type'].cat.categories

        # calculate centroids of cell types
        cellTypeCentroids = []
        for cellType in cellTypes:
            cellTypeCentroids.append(tuple(map(lambda y: sum(y) / float(len(y)), zip(*cellTypePointMap[cellType]))))

        # calculate SSE (sum of squared errors) across each cell type
        SSE_C = 0.0
        for cell_type_1 in range(len(cellTypes) - 1):
            for cell_type_2 in range(cell_type_1 + 1, len(cellTypes)):
                for coord in range(2):
                    SSE_C +=(cellTypeCentroids[cell_type_1][coord] - cellTypeCentroids[cell_type_2][coord]) ** 2

        # calculate MSE (mean squared error) across all cell types
        T = len(cellTypes)
        return SSE_C/(T*(T+1))

    def calculate_MSE_B(self):

        # create map of batches to points
        batchPointMap = {k: [] for k in self._adata.obs['batch'].cat.categories}
        for i in range(len(self._adata.obs)):
            batchPointMap[self._adata.obs.iloc[i]['batch']].append(self._adata.obsm['X_umap'][i])


        # create centroids of batches
        batchCentroids = {batch: tuple() for batch in self._adata.obs['batch'].cat.categories}
        for batch in self._adata.obs['batch'].cat.categories:
            batchCentroids[batch] = tuple(map(lambda y: sum(y) / float(len(y)), zip(*batchPointMap[batch])))

        # calculate SSE (sum of squared errors) for each batch
        SSE_B = {batch: 0.0 for batch in self._adata.obs['batch'].cat.categories}
        for batch in self._adata.obs['batch'].cat.categories:
            for point in batchPointMap[batch]:
                for coord in range(2):
                    SSE_B[batch] += (point[coord] - batchCentroids[batch][coord]) ** 2

        # calculate SSE (sum of squared errors) across all batches
        MSE_B = 0.0
        for batch in self._adata.obs['batch'].cat.categories:
            MSE_B += SSE_B[batch] / len(batchPointMap[batch])

        # calculate MSE (mean squared error) across all batches
        return (1/4) * MSE_B

    def calculate_F_statistic(self):

        # in this exercise, the F statistic is defined as the ratio of MSE_C to MSE_B
        return self.calculate_MSE_C() / self.calculate_MSE_B()


def main():
    '''
    this is an optional driver for the class
    '''
    if(len(sys.argv) != 2):
        sys.stderr.write("usage: " + __file__ + " <adata-file-path>\n")
        sys.exit(1)

    # read in adata object from file system
    adata = sc.read(sys.argv[1])

    # define number of neighbors k
    k = 12

    # calculate F statistic for batch-balanced

    # instantiate bblknn graph
    knnGraph = bblknn_graph(adata, k_per_batch=6, l=3, n_components=50)

    # run louvain to cluster data
    sc.tl.louvain(knnGraph._adata)

    # run umap to project in 2-space
    sc.tl.umap(knnGraph._adata)

    # instantiate and calculate F statistic
    myFStat = FStatistic(knnGraph._adata)
    print("batch-balanced f stat: " + str(myFStat.calculate_F_statistic()))


    # calculate F statistic for non-batch-balanced
    # read in AnnData object
    pbmc = sc.read(sys.argv[1])

    # instantiate knn graph
    myKNNG = KnnG(pbmc, n_neighbors=12)

    # run louvain to cluster data
    sc.tl.louvain(myKNNG._adata,  flavor='igraph', directed=False, use_weights=True)

    # reduce UMA
    sc.tl.umap(myKNNG._adata)

    # instantiate and calculate F statistic
    myFStat = FStatistic(myKNNG._adata)
    print("non-batch-balanced f stat: " + str(myFStat.calculate_F_statistic()))

if __name__ == "__main__":
    main()


