#! /usr/bin/env/ python

import numpy as np
import scipy
from euclid_knn import knnG
from sklearn.metrics import pairwise_distances
from scanpy.neighbors import compute_connectivities_umap
import logging

#EXAMPLE USAGE
# >>> clf = bbknn_graph(adata, batchlabel='label', neighbors_within_batch=6, pcs=50) #method is always umap
# >>> bbknn_indices, bbknn_distances = clf.bbknn() #this will compute batch balanced neighbors with your choice of neighbors_within_batch
# >>> l_k_bbknn_indices, l_k_bbknn_distances = clf.l_k_bbknn(l=3) # this will compute a subsample of your bbknn by using the knn indices and distances computed by clf.bbknn()


class bbknn_graph():
    '''
    Function performs batched balanced knn
    INPUT: 
        1. adata object with PCA done on it already
        2. batch label
        3. neighbors within each batch
    
    Output: new knn distances and connectivites with shape (n_observations, n_neighbors*n_batches)
    
    '''
    
    logger = logging.getLogger(__name__)

    def __init__(self,adata, batchlabel = None, neighbors_within_batch=6, 
                 pcs=50, method='umap', batch_unique=2):
        
        #fill in method
        self.batch_unique = batch_unique
        self.neighbors_within_batch = neighbors_within_batch
        
        #instantiating matrices for distances and indices
        if adata :
            self.knn_distances = np.zeros((adata.shape[0],neighbors_within_batch*len(self.batch_unique)))
            self.knn_indices = np.copy(self.knn_distances).astype(int)
        else:
            # unit test ubur
        
        # instantiating matrices for l-k-bbknn
            self.l_knn_indices = None
            self.l_knn_distances = None
        
        
        self.connectivities = None
        self.distances = None
    
    def get_connectivities(self,knn_indices, knn_distances):
        self.distances, self.connectivities = compute_connectivities_umap(knn_indices, knn_distances, knn_indices.shape[0], knn_indices.shape[1])
    
    def update_adata(self):
        #updating adata.uns
        self.adata.uns['neighbors']={}
        self.adata.uns['neighbors']['params'] = {}
        self.adata.uns['neighbors']['params']['n_neighbors']=self.neighbors_within_batch
        self.adata.uns['neighbors']['params']['method'] = self.method
    
        assert self.connectivities is not None
        assert self.distances is not None
        
        self.adata.uns['neighbors']['connectivities'] = self.connectivities
        self.adata.uns['neighbors']['distances'] = self.distances
    
    def querry(self,querry, ref):
        #return distances and indices of each new querry sample based on the knn obtained from get_knn
        # default is euclidean distance
    
        #fill in method
        pass
    
    def get_neighbors(self,D):
        ''' 
        function returns k most similar neighbors for each sample(row)
            Input is a distance matrix calculaed from the get_distances method
            Output is are two matrices:
                1. the distance of each sample (row) against every other sample (column) sorted from smallest distance (most similar) to greatest distance(least similar) for the k neighbors
                2. the index matrix of the sorted k nearest neighbors corresponding to the sorted distances above
        '''
    
        #fill in method
        pass
            
    
#     def bbknn(self):
#         '''
#         Function computes distances for all combinations of batches.
#         Fill in the below for loops
#         '''
#         for i in range(len(self.batch_unique)):
#             #get ref batch
#     
#     
#             for j in range(len(self.batch_unique)):
#     
#                 #querry new batch with ref batch
#                 #this means that the distances between the querry batch and the ref_batch are calculated
#     
#                 #Example:
#                 # if you have batches =3 and neighbors=3
#                 # your knn indices and distances matrices will have dims (pca_matrix.shape[0], n_batches*neighbors) = (n_obs, 9)
#                 # in order ot update these matrices you need to do the following:
#                 # for the first batch update the 0-k (n_neighbbors) columns
#                 # for the second batch update the k-2*k (k=n_neighbors) columns
#                 # for the third batch update the 2*k - 3*k (k=n_neighbors) columns
#                 # the indeces for the rows are always what you're querrying
#                 
#                 #the first k columns are batch1-batch1 and batch1-batch2
#                 #the next k columns are batch2-batch1 and batch2-batch2
                
    def l_k_bbknn(self,l=2):
        #this method makes an l subsampling of the bbknn computed in the graph() method
        # if k=4, meaning you are finding 4 neighbors between batches, and l=2:
                ## subsample 2 neighbors from batch1 and 2 neighbors from batch2
            
            
        if l >= self.neighbors_within_batch:
            raise ValueError('l cannot be equal or larger than k')
            
        self.l_knn_indices = np.zeros((self.knn_indices.shape[0], 2*l)).astype(int)
        self.l_knn_distances = np.zeros((self.knn_distances.shape[0], 2*l))
            
        nCols = self.knn_indices.shape[0]
        self.logger.info("nCols:{}".format(nCols))
        
        #         for i in range(len(self.batch_unique)):
        #             #fill in loop
        #             pass
        for rowIdx in range(self.l_knn_indices.shape[0]):
            # select random index for each k groups
            # if k=3, self.knn_indices has 6 columns
            # so here we are going to select l random indices from first three columns, and second three columns
            # for the above example i am assuming two neighbors
                
            # split into batch_unique arrays
            rowIndices = self.knn_indices[rowIdx,:]
            rowDist = self.knn_distances[rowIdx, :]
            
            self.logger.info("rowIndices:{}".format(rowIndices))
            self.logger.info("rowDist   :{}".format(rowDist))

            
            rowIndicesSplits = np.split(rowIndices, self.batch_unique)
            rowDistSplits = np.split(rowDist, self.batch_unique)
            
            self.logger.info("rowIndicesSplits:\n{}".format(rowIndicesSplits))
            self.logger.info("rowDistSplits:\n{}".format(rowDistSplits))
            
            tmpIndices = np.zeros(self.batch_unique * l, dtype=int)
            tmpDistances = np.zeros(self.batch_unique * l)
            
            self.logger.info("***********8")
            for splitIdx in range(len(rowIndicesSplits)):
                low = 0 # splitIdx * self.batch_unique
                high = self.batch_unique #, low + self.batch_unique
                rIdx = np.random.randint(low, high, l)
                
                self.logger.info("l:{} low:{} high:{}  rIdx:{}".format(l, low, high, rIdx))
                self.logger.info("splitIdx:{}".format(splitIdx))
                
                start = splitIdx * l
                end = start + l
                self.logger.info("start:{} end:{}\n".format(start, end))
                tmpIndices[start:end]   = rowIndicesSplits[splitIdx][rIdx]
                tmpDistances[start:end] = rowDistSplits[splitIdx][rIdx]
                
            
            self.l_knn_indices[rowIdx]   = tmpIndices
            self.l_knn_distances[rowIdx] = tmpDistances
            
            self.logger.info("%%%% tmpIdx:{}".format(tmpIndices))
            self.logger.info("%%% tmpIdx:{}\n".format(tmpDistances))






