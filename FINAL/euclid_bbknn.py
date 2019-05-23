#! /usr/bin/env/ python

from euclid_knn import KnnG
import logging
import numpy as np
# import scipy
from scanpy.neighbors import compute_connectivities_umap
import sys
import scanpy as sc


# from scipy import sparse

################################################################################
class bbknn_graph():
    '''
    Output: new knn distances and connectivites with shape (n_observations, n_neighbors*n_batches)

    public functions:
        def __init__(self,adata, batchLabel = None, neighbors_within_batch=6,
                 runPCA =True, pcs=50, method='umap', batch_unique=2)

        def l_k_bbknn(self,l=2):

    updates:
        adata.obsm['X_pca']
        adata.uns['neighbors']['connectivities']
        adata.uns['neighbors']['distances']
    '''

    logger = logging.getLogger(__name__)

    ######################################################################
    def __init__(self, adata, batchLabel=None, neighbors_within_batch=6,
                 runPCA=True, pcs=50, method='umap'):
        '''
        input:

        '''

        self._adata = adata
        self._batchlabel = batchLabel
        self._neighbors_within_batch = neighbors_within_batch
        self._runPCA = runPCA
        self._psc = pcs
        self._method = method
        self._numBatches = None

        self._knn_distances = None
        self._knn_indices = None
        self._l_knn_indices = None
        self._l_knn_distances = None
        self._connectivities = None
        self._distances = None

        if adata:
            batchCounts = self._calcNumCellsInEachBatch()
            self._numBatches = len(batchCounts)
            numRows = adata.shape[0]
            numCols = neighbors_within_batch * self._numBatches
            self._knn_distances = np.zeros((numRows, numCols))
            self._knn_indices = np.zeros((numRows, numCols), dtype=int)

            # run KnnG it will run pca and calculate nei
            D = self._calcPairWiseDistanceMatrix()

            self._knn_indices, self._knn_distances = self._bbknn(D, batchCounts)
            self._get_connectivities(self._knn_indices, self._knn_distances)
            self._update_adata()
        else:
            # unit test
            pass

    ######################################################################
    def _bbknn(self, D, batchCounts):
        '''
        use Knng to compute nearest Neighbors for each batch and combine the results

        returns balanced knn_indices,knn_dist

        input:
            D: pair wise distance matrix

            batchCounts:
                pass none in production
                unit test should pass something like [('0', 8098), ('1', 7378)]
        '''
        knng = KnnG(None)  # init like unit test to reduce run time

        knng._n_neighbors = self._neighbors_within_batch

        # split D up by batches
        if not batchCounts:
            batchCounts = self._calcNumCellsInEachBatch()

        splitsLocations = self._calcSplits(batchCounts)

        #
        # we split by cols
        # explanation: assume we have two batches with different number of cells
        # when we calculate the nearest neighbors using the first split
        # the results will be [[batch1 x batch 1], [batch2 x batch1]]
        #
        # when we use the second batch we get
        # [[batch1 x batch2], [batch2 x batch2]]
        #
        byCols = 1
        splits = np.split(D, splitsLocations, axis=byCols)
        # np.split(D, [3,6] returns
        # [:3], [3:6], [6:]
        # we need to remove this last split. It is empty
        del splits[-1]

        # for each batch calculate the nearest neighbors
        batchNNIdx = []
        batchNNDist = []
        for i in range(len(splits)):
            split = splits[i]
            batchIdx, batchDist = knng._get_neighbors(split)

            # we need to adj the indexes so that they map back to the columns
            # in our original matrix
            offset = self._calcBatcholOffset(splitsLocations, i)
            batchIdx = batchIdx + offset
            self.logger.info("batchIdx.shape:{} batchDist.shape{}".format(batchIdx.shape, batchDist.shape))

            #             self.logger.info("batchIdx:\n{}".format(batchIdx))
            #             self.logger.info("batchNNDist:\n{}".format(batchNNDist))
            batchNNIdx.append(batchIdx)
            batchNNDist.append(batchDist)

        # concatenate the batches to create the balanced batch nearest neighbors
        byCols = 1
        knn_indices = np.concatenate(batchNNIdx, axis=byCols)
        knn_dist = np.concatenate(batchNNDist, axis=byCols)

        return knn_indices, knn_dist

    ######################################################################
    def _calcBatcholOffset(self, splitsLocations, i):
        '''
        np.split() returns matrix slices. the indexing starts at zero
        we need to calculate the offset back to the original matrix
        '''
        ret = 0
        if i > 0:
            ret = splitsLocations[i - 1]

        return ret

    ######################################################################
    def _calcNumCellsInEachBatch(self):
        '''
        returns
            example [('0', 8098), ('1', 7378)]
        '''
        df = self._adata.obs[['batch']]
        batchIdx = df['batch'].unique()
        batchKeys = [i for i in batchIdx]

        ret = []
        for bk in batchKeys:
            # print("bk:{} type(bk):{}".format(bk, type(bk)))
            rows = df.loc[:, ['batch']] == bk  # int(bk)
            # print(sum(rows.values))
            ret.append((bk, sum(rows.values)[0]))

        return ret

    ######################################################################
    def _calcPairWiseDistanceMatrix(self):
        '''
        use KnnG to compute the pair wise distance
        we are in the same package and make use of non public parts of the implement
        '''
        knng = KnnG(None)  # init like a unit test to reduce run time
        knng._adata = self._adata
        knng._PCA(npc=50)
        knng._calDistance()

        return knng._D

    ######################################################################
    def _calcSplits(self, batchCounts):
        '''
        calculates how to split D base on batch sizes

        return example [3, 7]
        '''
        splits = []
        start = 0
        for i in range(len(batchCounts)):
            bk, bc = batchCounts[i]
            splits.append(start + bc)
            start += bc

        self.logger.info("splits:{}".format(splits))

        return splits

    ######################################################################
    def _get_connectivities(self, knn_indices, knn_distances):
        d, c = compute_connectivities_umap(knn_indices,
                                           knn_distances,
                                           knn_indices.shape[0],
                                           knn_indices.shape[1])
        self._distances = d
        self._connectivities = c


    ######################################################################
    def _update_adata(self, l=None):
        '''
        set up data correctly for scanpy.louvain:
        if l_k_bbnn
            'n_neighbors' = self._numBatches * l
        else bbknn
            'n_neighbors' = self._neighbors_within_batch
        '''

        # updating adata.uns
        self._adata.uns['neighbors'] = {}
        self._adata.uns['neighbors']['params'] = {}

        n = self._neighbors_within_batch
        if l:
            n = self._numBatches * l  # (l - 1) aedwip do not subtrack -1
            self.logger.info('n:{} l :{}, self._numBatches: {}'.format(n, l, self._numBatches))

        self._adata.uns['neighbors']['params']['n_neighbors'] = n
        self._adata.uns['neighbors']['params']['method'] = self._method

        assert self._connectivities is not None
        assert self._distances is not None

        self._adata.uns['neighbors']['connectivities'] = self._connectivities
        self._adata.uns['neighbors']['distances'] = self._distances

def main():
    '''
    this is an optional driver for the class
    '''
    if(len(sys.argv) != 2):
        sys.stderr.write("usage: " + __file__ + " <adata-file-path>\n")
        sys.exit(1)

    # read in adata object from file system
    adata = sc.read(sys.argv[1])

    # build bblknn graph
    myGraph = bbknn_graph(adata=adata, batchLabel=None, neighbors_within_batch=6, runPCA=True, pcs=50, method='umap')

    # run louvain to cluster data
    sc.tl.louvain(myGraph._adata)

    # run umap to project in 2-space
    sc.tl.umap(myGraph._adata)

    # plot the bbknn graph
    sc.pl.umap(myGraph._adata,  color=['Cell type', 'batch', 'louvain'])

if __name__ == "__main__":
    main()