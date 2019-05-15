#! /usr/bin/env/ python

from euclid_knn import KnnG
import logging
import numpy as np
import pandas as pd
#import scipy
from sklearn.metrics import pairwise_distances
from scanpy.neighbors import compute_connectivities_umap


################################################################################
class bbknn_graph():
    '''
    TODO: AEDWIP:
    
    Function performs batched balanced knn
    INPUT: 
        1. adata object with PCA done on it already
        2. batch label
        3. neighbors within each batch
    
    Output: new knn distances and connectivites with shape (n_observations, n_neighbors*n_batches)
    
    public functions:
        def __init__(self,adata, batchLabel = None, neighbors_within_batch=6, 
                 runPCA =True, pcs=50, method='umap', batch_unique=2)
                 
        def l_k_bbknn(self,l=2):
    '''
    
    logger = logging.getLogger(__name__)
    
    ######################################################################
    def __init__(self,adata, batchLabel = None, neighbors_within_batch=6, 
                 runPCA =True, pcs=50, method='umap', batch_unique=2):
        '''
        input:
            batch_unique: 
                the number of batches in the combined data set adata.X.
                assume the batches are stack in continuously blocks horizontally
            
        '''
        
        self._adata = adata
        self._batchlabel = batchLabel
        self._neighbors_within_batch = neighbors_within_batch
        self._runPCA = runPCA
        self._psc = pcs
        self._method = method
        self._batch_unique = batch_unique 
        
        self._knn_distances = None
        self._knn_indices = None
        self.l_knn_indices = None
        self.l_knn_distances = None
        self._connectivities = None
        self._distances = None
        
        if adata :
            self._knn_distances = np.zeros((adata.shape[0],neighbors_within_batch*len(self._batch_unique)))
            self._knn_indices = np.copy(self.knn_distances).astype(int)
            
            # run KnnG it will run pca and calculate nei
            D = self._calcPairWiseDistanceMatrix()
            
            knn_indices, knn_dist = self._bbknn(D)
            self._get_umap_connectivities(knn_indices, knn_dist)
            self._update_adata()
        else:
            # unit test 
            pass

    ######################################################################
    def _bbknn(self, D):
        '''
        use Knng to compute nearest Neighbors for each batch and combine the results
        
        returns balanced knn_indices,knn_dist
        '''
        knng = KnnG(None) #init like unit test to reduce run time
        
        knng._n_neighbors = self._neighbors_within_batch * self._batch_unique
        
        # split D up by batches
        batchCounts = self._calcNumCellsInEachBatch()
        splitsLocations = self._calcSplits(batchCounts)
        #
        # we split by cols
        # explanation: assume we have two batchs with different number of cells
        # when we calculate the nearest neighbors using the first split
        # the results will be [[batch1 x batch 1], [batch2 x batch1]]
        #
        # when we use the second batch we get
        # [[batch1 x batch2], [batch2 x batch2]]
        #
        byCols = 1
        splits = np.split(D, splitsLocations, axis=byCols)
        
        # for each batch calculate the nearest neighbors
        batchNNIdx = []
        batchNNDist = []
        for split in splits:
            self.logger.info("split.shape:{}".format(split.shape))
            batchIdx, batchDist =  knng._get_neighbors(split)
            self.logger.info("batchIdx.shape:{} batchDist.shape{}".format(batchIdx.shape, batchDist.shape))
            batchNNIdx.append(batchIdx)
            batchNNDist.append(batchDist)
            
        # concatenate the batches to create the balanced batch nearest neighbors
        byRows = 0
        knn_indices  = np.concatenate(batchNNIdx, axis=byRows)
        knn_dist = np.concatenate(batchNNDist, axis=byRows)
        
        return knn_indices,knn_dist

    ######################################################################                
    def _calcSplits(self, batchCounts):
        '''
        calculates how to split D base on batch sizes
        '''
        splits = []
        start = 0
        for i in range(len(batchCounts)):
            bk, bc = batchCounts[i]
            splits.append(start + bc)
            start += bc
        self.logger.info("splits:{}".format(splits))        
        
        
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
            #print("bk:{} type(bk):{}".format(bk, type(bk)))
            rows = df.loc[:, ['batch']] == bk #int(bk)
            print(sum(rows.values))
            ret.append((bk, sum(rows.values)[0]))
            
        return ret
    
    ######################################################################                
    def _calcPairWiseDistanceMatrix(self):
        '''
        use KnnG to compute the pair wise distance 
        we are in the same package and make use of non public parts of the implement
        '''
        knng = KnnG(None) # init like a unit test to reduce run time
        knng._adata = self._adata
        knng._PCA(npc=50)
        knng._calDistance()
        
        return knng._D
        
    ######################################################################    
    def _get_connectivities(self,knn_indices, knn_distances):
        d,c = compute_connectivities_umap(knn_indices, 
                                          knn_distances, 
                                          knn_indices.shape[0], 
                                          knn_indices.shape[1])
        self._distances = d
        self._connectivities = c
    
    ######################################################################    
    def _update_adata(self):
        #updating adata.uns
        self._adata.uns['neighbors']={}
        self._adata.uns['neighbors']['params'] = {}
        self._adata.uns['neighbors']['params']['n_neighbors']=self._neighbors_within_batch
        self._adata.uns['neighbors']['params']['method'] = self._method
    
        assert self.connectivities is not None
        assert self.distances is not None
        
        self._adata.uns['neighbors']['connectivities'] = self._connectivities
        self._adata.uns['neighbors']['distances'] = self._distances
    
#     ######################################################################    
#     def querry(self,querry, ref):
#         #return distances and indices of each new querry sample based on the knn obtained from get_knn
#         # default is euclidean distance
#     
#         #fill in method
#         pass
    
#     ######################################################################
#     def get_neighbors(self,D):
#         ''' 
#         function returns k most similar neighbors for each sample(row)
#             Input is a distance matrix calculaed from the get_distances method
#             Output is are two matrices:
#                 1. the distance of each sample (row) against every other sample (column) sorted from smallest distance (most similar) to greatest distance(least similar) for the k neighbors
#                 2. the index matrix of the sorted k nearest neighbors corresponding to the sorted distances above
#         '''
#     
#         #fill in method
#         pass
            
    
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
           
    ######################################################################
    def l_k_bbknn(self,l=2):
        '''
        this method makes an l subsampling of the bbknn computed in the graph() 
        method. if k=4, meaning you are finding 4 neighbors between batches, 
        and l=2: subsample 2 neighbors from batch1 and 2 neighbors from batch2
        '''
        if l >= self._neighbors_within_batch:
            raise ValueError('l cannot be equal or larger than k')
            
        # init storage for results
        self.l_knn_indices = np.zeros((self.knn_indices.shape[0], 2*l)).astype(int)
        self.l_knn_distances = np.zeros((self.knn_distances.shape[0], 2*l))
            
        nCols = self.knn_indices.shape[0]
        self.logger.debug("nCols:{}".format(nCols))
        
        for rowIdx in range(self.l_knn_indices.shape[0]):
            # select the current rows to process   
            rowIndices = self.knn_indices[rowIdx,:]
            rowDist = self.knn_distances[rowIdx, :]
            
            # split into batch_unique number of arrays   
            rowIndicesSplits = np.split(rowIndices, self._batch_unique)
            rowDistSplits = np.split(rowDist, self._batch_unique)
                        
            # init storage for tmp results
            tmpIndices = np.zeros(self._batch_unique * l, dtype=int)
            tmpDistances = np.zeros(self._batch_unique * l)
            
            for splitIdx in range(len(rowIndicesSplits)):
                # randomly select index values 
                low = 0 
                high = self._batch_unique 
                rIdx = np.random.randint(low, high, l)
                            
                # copy split selections to final results    
                start = splitIdx * l
                end = start + l
                tmpIndices[start:end]   = rowIndicesSplits[splitIdx][rIdx]
                tmpDistances[start:end] = rowDistSplits[splitIdx][rIdx]
                
            
            self.l_knn_indices[rowIdx]   = tmpIndices
            self.l_knn_distances[rowIdx] = tmpDistances
            
            self.logger.debug("tmpIdx:{}".format(tmpIndices))
            self.logger.debug("tmpIdx:{}\n".format(tmpDistances))
            
        # TODO: AEDWIP:
        self.logger.error("AEDWIP: need to call _get_connectivities() and _update()")

