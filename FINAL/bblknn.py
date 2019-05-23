from sklearn.metrics import pairwise_distances
from sklearn.decomposition import PCA
import random
import scanpy as sc
import numpy as np
from euclid_bbknn import bbknn_graph
import logging
import sys


################################################################################
class bblknn_graph():
    '''
    Output: new knn distances and connectivites with shape (n_observations, n_neighbors*n_batches)

    methods:
        def __init__(self,adata, k_per_batch, l, n_components)

        def l_k_bbknn(self,l=2):

    updates:
        adata.obsm['X_pca']
        adata.uns['neighbors']['connectivities']
        adata.uns['neighbors']['distances']
    '''

    logger = logging.getLogger(__name__)

    ######################################################################
    def __init__(self, adata, k_per_batch, l, n_components):
        '''
        input:
            adata (AnnData object)

            k_per_batch (number of neighbors per batch)

            l (number of neighbors to subsample from each batch

            n_components (number of principle components)

        '''

        # validate values of l and k
        if l >= k_per_batch:
            raise ValueError('l must be strictly less than k')

        # initialize member variables
        self._adata = adata
        self._neighbors_within_batch = k_per_batch
        self._l = l
        self._n_components = n_components

        # run l_k_bbknn on adata
        self.l_k_bbknn()

    def _runPCA(self):
        '''
        this method is a simple wrapper for sklearn PCA method
        '''

        # find principle components of adata using PCA
        pca = PCA(self._n_components)

        # set adata X_pca data structure to transformed adata
        self._adata.obsm['X_pca'] = pca.fit_transform(self._adata.X)

    def l_k_bbknn(self):

        '''
        this method makes an l subsampling of the batch-balanced neighbors from bbknn

        if k=6 and l=3, then subsample 3 neighbors from batch1 and 3 neighbors from batch2
        '''

        # get list of batch identifiers
        batch_unique = self._adata.obs.batch.cat.categories

        # run pca
        self._runPCA()

        # create bbknn_graph
        bbknn = bbknn_graph(self._adata, neighbors_within_batch= self._neighbors_within_batch, runPCA=False, pcs=self._n_components)

        # create indices for random sampling l from k
        random_sample = random.sample(range(0, self._neighbors_within_batch), self._l)

        # initialize lk_bbknn l_bbknn_indices and l_bbknn_distances matrices to 0's
        l_bbknn_indices = np.zeros((bbknn._knn_indices.shape[0], len(batch_unique) * self._l)).astype(int)
        l_bbknn_distances = np.zeros((bbknn._knn_distances.shape[0], len(batch_unique) * self._l))

        # outer loop through batches
        for i in range(len(batch_unique)):
            # get batch id for ref_batch
            batch_id = batch_unique[i]

            # get booleen index for reference batch
            bool_idx = self._adata.obs['batch'] == batch_id

            # use booleen index to get pca data for reference batch
            ref_batch_pca = self._adata.obsm['X_pca'][bool_idx]

            # create a booleen index for ref_batch to map back to pca matrix
            ref_batch_idx = np.arange(self._adata.shape[0])[bool_idx]

            # inner loop through batches
            for j in range(len(batch_unique)):
                # get batch id for query_batch
                batch_id = batch_unique[j]

                # get booleen index for query batch
                bool_idx = self._adata.obs['batch'] == batch_id

                # use booleen index to get pca data for query batch
                query_batch_pca = self._adata.obsm['X_pca'][bool_idx]

                # create a booleen index for query_batch to map back to pca matrix
                query_batch_idx = np.arange(self._adata.shape[0])[bool_idx]

                # calculate pairwise_distances between query batch and ref_batch
                D = pairwise_distances(X=query_batch_pca, Y=ref_batch_pca)

                # get indices for n nearest neighbors
                neighbors = np.argsort(D, axis=1)[0:,0:self._neighbors_within_batch]

                # get distance for n nearest neighbors (including self)
                sorted_D = np.sort(D, axis=1)[0:,0:self._neighbors_within_batch]

                # map nearest neighbors to pca indices
                for n in range(neighbors.shape[0]):
                    for k in range(neighbors.shape[1]):
                        temp_neighbor = neighbors[n,k]
                        neighbors[n,k]=ref_batch_idx[temp_neighbor]

                # set range of columns for indices and distances
                col_range = np.arange(i * self._l, (i+1) * self._l)

                # pass random sampled l nearest neighbors to l_bbknn_indices and distances matrix
                l_bbknn_indices[query_batch_idx[:,None], col_range[None,:]]= neighbors[:,random_sample]
                l_bbknn_distances[query_batch_idx[:,None], col_range[None,:]] = sorted_D[:,random_sample]

        # calculate connectivities and distances using scanpy method
        distances, connectivities = sc.neighbors.compute_connectivities_umap(l_bbknn_indices, l_bbknn_distances, n_obs=len(self._adata.obs),n_neighbors=2 * self._l)

        # set connectivities and distances in adata object
        self._adata.uns['neighbors']['connectivities'] = connectivities
        self._adata.uns['neighbors']['distances'] = distances


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
    myGraph = bblknn_graph(adata, k_per_batch=6, l=3, n_components=50)

    # run louvain to cluster data
    sc.tl.louvain(myGraph._adata)

    # run umap to project in 2-space
    sc.tl.umap(myGraph._adata)

    # plot the bblknn graph
    sc.pl.umap(myGraph._adata,  color=['Cell type', 'batch','louvain'])

if __name__ == "__main__":
    main()

