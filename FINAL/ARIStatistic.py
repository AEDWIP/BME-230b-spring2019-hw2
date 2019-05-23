from sklearn.metrics.cluster import adjusted_rand_score
import logging
from bblknn import bblknn_graph
import scanpy as sc
from matplotlib import pyplot
import numpy as np
from euclid_bbknn import bbknn_graph
import sys




################################################################################
class ARIStatistic():
    '''
    Output: list of adjusted rand indices for n_samples samples

    methods:
        __init__(self,adata, k_per_batch, l, n_components)

        calculate_ARI(): calculate adjusted rand index
    '''

    logger = logging.getLogger(__name__)

    ######################################################################
    def __init__(self, adata, k_per_batch, l, n_components, n_samples, bbknn_louvain):
        '''
        input:
            adata (AnnData object)

            k_per_batch (number of neighbors per batch)

            l (number of neighbors to subsample from each batch

            n_components (number of principle components)

            n_samples (number of samples to run)

            bbknn_louvain (louvain clustering of batch-balanced knn data)

        '''

        # number of samples to run
        self._samples = n_samples

        # AnnData object
        self._adata = adata

        # number of neighbors per batch
        self._k_per_batch = k_per_batch

        # number of neighbors to subsample
        self._l = l

        # number of principle components
        self._n_components = n_components

        # louvain clustering of batch-balanced knn data
        self._bbknn_louvain = bbknn_louvain

        # list of ARIs obtained from each subsample
        self._results = []

        for sample in range(self._samples):
            self._adata.obs['louvain'] = None

            mybblknn = bblknn_graph(adata, k_per_batch=self._k_per_batch, l=self._l, n_components=self._n_components)

            # run louvain to cluster data
            sc.tl.louvain(mybblknn._adata)

            # run umap to project in 2-space
            sc.tl.umap(mybblknn._adata)

            # run louvain clustering
            sc.tl.louvain(mybblknn._adata)

            # calculate adjust rand index for sample
            ars = adjusted_rand_score(mybblknn._adata.obs['louvain'], self._bbknn_louvain)

            # append ari to results list
            self._results.append(ars)
            
            print("finished sample " + str(sample))
            print("ars = " + str(ars))


def main():
    '''
    this is an optional driver for the class
    '''
    if(len(sys.argv) != 2):
        sys.stderr.write("usage: " + __file__ + " <adata-file-path>\n")
        sys.exit(1)

    # read in adata object from file system
    adata_bbknn = sc.read(sys.argv[1])

    # run bbknn on adata
    bbknn = bbknn_graph(adata_bbknn, neighbors_within_batch=6, runPCA=True, pcs=50)

    # run louvain to cluster data
    sc.tl.louvain(bbknn._adata)

    # read in AnnData object (again)
    adata_blknn = sc.read(sys.argv[1])

    # instantiate ARIStatistic object
    myARIStat = ARIStatistic(adata=adata_blknn, k_per_batch=6, l=3, n_components=50, n_samples=10,
                             bbknn_louvain=bbknn._adata.obs['louvain'])

    # plot ARIs
    pyplot.bar(range(1, len(myARIStat._results) + 1), myARIStat._results, color='black')

    # print statistics
    print('Average:', np.mean(myARIStat._results), 'SD:', np.std(myARIStat._results))

if __name__ == "__main__":
    main()