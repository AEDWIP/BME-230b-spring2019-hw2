#
# BME-230B Spring 2019 HW 2 Question 1
# James Casaletto, Andrew Davidson, Yuanqing Xue, Jim Zheng
#
#
# Question 1.a, 1.b see [euclid_knn.py](euclid_knn.py)
#
#
# ref:
# - [ Single-Cell Analysis in Python](https://scanpy.readthedocs.io/en/stable/api/index.html#tools-tl)
# - [data exploration](exploreData.ipynb)
# - [scanpy.tl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.umap.html)
# - <span style="color:red">scanpy.api.pl no longer exists</span>
# - [scanpy.pl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.pl.umap.html#scanpy.pl.umap)

import logging
import numpy as np
from scanpy.neighbors import compute_connectivities_umap
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
import sys
import scanpy as sc


################################################################################
class KnnG():
    '''
    replaces scanpy pca and k nearest  neighbors  with our own

    updates:
        adata.obsm['X_pca']
        adata.uns['neighbors']['connectivities']
        adata.uns['neighbors']['distances']

    usage:
        knn = KnnG(adata)
        scanpy.tl.umap(adata)
        scanpy.pl.umap(adata)

    public functions
        __init__(self, adata = None, d_metric='euclidean', n_neighbors=15, method='umap',
                    runPCA=True, nPC=50)
    '''
    logger = logging.getLogger(__name__)

    ######################################################################
    def __init__(self, adata=None, d_metric='euclidean', n_neighbors=15,
                 method='umap', runPCA=True, nPC=50):
        '''
        arguments:
            adata: unit test should pass None and init members as needed
            nPC: if runPCA == True should be the number of principle components you want
                to reduce adata.X to
        '''
        self._adata = adata
        self._d_metric = d_metric
        self._n_neighbors = n_neighbors
        self._method = method

        self._D = None  # the pair wise distance a numpy matrix.
        self._connectivities = None
        self._distances = None

        # wrapped initialization to make unit test run faster an easier to write
        # unit test should construct KnnG(adata=None)
        if self._adata:
            self.logger.info("calculating please be patient")
            self._adata.uns['neighbors']['connectivities'] = None
            self._adata.uns['neighbors']['distances'] = None

            # calculate k neighbors and umap connectivities:
            if runPCA:
                self._PCA(nPC)

            self._calDistance()
            knn_indices, knn_dist = self._get_neighbors(self._D)
            self._get_umap_connectivities(knn_indices, knn_dist)

            self._update_adata()

    ######################################################################
    def _calDistance(self):
        '''
        calculates pair wise distances. results store in self.D

        input:
        '''
        self.logger.info("BEGIN")
        self._D = pairwise_distances(self._adata.obsm['X_pca'])
        self.logger.info("END\n")

    ######################################################################
    def _findNeigborsForRow(self, row, k):
        '''
        arguments:
            row: a row in the pair wise distance matrix
            k: the number of nearest neighbors to find

        returns:
            a row in adjacency matrix format of k nearest neighbors as
            two numpy arrays. The values in the first array are the indices
            for the row argument nearest neighbors. The values of the second
            array are the distances
        '''

        # create fast way to sort distances and look up
        # corresponding distances
        distanceReverseIndex = {row[i]: i for i in range(len(row))}

        distances = distanceReverseIndex.keys()

        #
        # if we are running using results for true knn the diagonal will be 0
        # if we are running bb-knn some the distance matrix is not square
        # that is to say we may or may not have a zero
        #
        sortedDistances = sorted(distances)
        start = 0
        if sortedDistances[0] == 0:
            start = 1

        neighborsDistances = sorted(distances)[start: k + start]

        retIdx = np.zeros(k)
        retDist = np.zeros(k)
        for i in range(len(neighborsDistances)):
            distance = neighborsDistances[i]
            idx = distanceReverseIndex[distance]
            retIdx[i] = idx
            retDist[i] = distance

        return retIdx, retDist

    ######################################################################
    def _get_neighbors(self, D):
        '''
        arguments
            D: pairwise distance matrix

        returns:
            a row in adjacency matrix format of k nearest neighbors as
            two numpy arrays. The values in the first array are the indices
            for the row argument nearest neighbors. The values of the second
            array are the distances
        '''
        self.logger.info("BEGIN")
        n = D.shape[0]
        knn_i = np.zeros((n, self._n_neighbors))
        knn_d = np.zeros((n, self._n_neighbors))
        for i in range(n):
            row = D[i, :]
            neigborsIdx, neigborsDist = self._findNeigborsForRow(row, self._n_neighbors)
            knn_i[i] = neigborsIdx
            knn_d[i] = neigborsDist

        self.logger.info("END\n")
        return knn_i, knn_d

    ######################################################################
    def _get_umap_connectivities(self, knn_indices, knn_dists):
        '''
        ref:
            https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.Neighbors.html?highlight=scanpy.neighbors

        results are store in
            self.distances = distances
            self.connectivities = connectivities
        '''
        self.logger.info("BEGIN")
        # you have to read the code to figure out how to call compute_connectivities_umap
        n_obs = self._adata.X.shape[0]

        self.logger.info("type(knn_dists):{} knn_dists.shape:{}".format(type(knn_dists), knn_dists.shape))
        self.logger.info("knn_dists[0:3:\n{}".format(knn_dists[0:3]))

        self.logger.info("type(knn_indices):{}".format(type(knn_indices)))
        # does not have 'shape' self.logger.info("knn_indices.shape:{}".format(knn_indices.shape)
        self.logger.info("knn_indices[0:3]:\n{}".format(knn_indices[0:3]))
        # ???

        # knn_indices is ?? array of rows. each contains the idx of its nearest neighbors
        # knn_dists is a similar structure, the values are the distance to the neighboors
        distances, connectivities = compute_connectivities_umap(knn_indices,
                                                                knn_dists,
                                                                n_obs,
                                                                self._n_neighbors)

        self.logger.info("type(distances):{}".format(type(distances)))
        self.logger.info("distances:\n{}".format(distances))

        self.logger.info("type(connectivities):{}".format(type(connectivities)))
        self.logger.info("connectivities:\n{}".format(connectivities))

        self.logger.info("END")

        self._distances = distances
        self._connectivities = connectivities

    ######################################################################
    def _PCA(self, npc=15):
        '''
        input:
            npc: the n_components argument. should be either
                0 <= n_components <= 0
                or
                an int greater than 0

                if <= 0, n_components specifies the amount of variation to preserve
                else it is the number of dimension to reduce to

        results are stored in self.adata.obsm['X_pca']
        ref:
            - https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
        '''
        self.logger.info("BEGIN input shape:{} pc:{}".format(self._adata.shape, npc))

        pcaTransformer = PCA(n_components=npc)
        self._adata.obsm['X_pca'] = pcaTransformer.fit_transform(self._adata.X)

        self.logger.info("END return shape:{}".format(self._adata.obsm['X_pca'].shape))

        ######################################################################

    def _update_adata(self):
        '''
        The trick to replacing scanpy implementation with our own is to
        update the anndata object with our intermediate values

        we replace the scanpy version of PCA by updating
        adata.obsm['X_pca'] = our PCA() output

        we replace the scanpy version of k-nearest-neighbors by updating
        self.adata.uns['neighbors']['connectivities'] = our knn() output
        self.adata.uns['neighbors']['distances'] = out knn() output
        '''
        self.logger.info('BEGIN')
        self._adata.uns['neighbors'] = {}
        self._adata.uns['neighbors']['params'] = {}
        self._adata.uns['neighbors']['params']['n_neighbors'] = self._n_neighbors
        self._adata.uns['neighbors']['params']['method'] = self._method

        self._adata.uns['neighbors']['connectivities'] = self._connectivities
        self._adata.uns['neighbors']['distances'] = self._distances

        self.logger.info('END\n')


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
    myGraph = KnnG(adata=adata, d_metric='euclidean', n_neighbors=15, method='umap', runPCA=True, nPC=50)

    # run louvain to cluster data
    sc.tl.louvain(myGraph._adata)

    # run umap to project in 2-space
    sc.tl.umap(myGraph._adata)

    # plot the knn graph
    sc.pl.umap(myGraph._adata,  color=['Cell type', 'batch', 'louvain'])

if __name__ == "__main__":
    main()
