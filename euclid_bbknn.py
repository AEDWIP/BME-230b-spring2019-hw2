#! /usr/bin/env/ python

import numpy as np
import scipy
from euclid_knn import knnG
from sklearn.metrics import pairwise_distances
from scanpy.neighbors import compute_connectivities_umap
from scanpy.tools._utils import doc_n_pcs
from nbformat.v4.tests.nbexamples import cells
from scipy.spatial.distance import pdist, squareform

#EXAMPLE USAGE
# >>> clf = bbknn_graph(adata, batchlabel='label', neighbors_within_batch=6, pcs=50) #method is always umap
# >>> bbknn_indices, bbknn_distances = clf.bbknn() #this will compute batch balanced neighbors with 
#            your choice of neighbors_within_batch
# >>> l_k_bbknn_indices, l_k_bbknn_distances = clf.l_k_bbknn(l=3) # this will compute a 
#            subsample of your bbknn by using the knn indices and distances computed by clf.bbknn()


class bbknn_graph():
    '''
    Function performs batched balanced knn
    INPUT: 
        1. adata object with PCA done on it already
        2. batch label
        3. neighbors within each batch
    
    Output: new knn distances and connectivites with shape (n_observations, n_neighbors*n_batches)
    
    ref:
    - [Fast Batch Alignment of Single Cell
Transcriptomes Unifies Multiple Mouse Cell
Atlases into an Integrated Landscape](https://www.biorxiv.org/content/biorxiv/early/2018/08/22/397042.full.pdf)
    
    '''
    def __init__(self,adata, batchlabel = None, neighbors_within_batch=6, pcs=50, method='umap'):
        #fill in method
        self.adata = adata
        self.batchlabel = batchlabel
        self.neighbors_within_batch = neighbors_within_batch
        self.pcs=doc_n_pcs
        self.method = method
        
        #instantiating matrices for distances and indices
        self.knn_distances = np.zeros((adata.shape[0],neighbors_within_batch*len(self.batch_unique)))
        self.knn_indices = np.copy(self.knn_distances).astype(int)
        #instantiating matrices for l-k-bbknn
        self.l_knn_indices = None
        self.l_knn_distances = None
        
        # compute pairwise distance
        # see doc in euclid_knn get_distance()
        condensedDistances = pdist(self.adata.obsm['X_pca'], metric=self.d_metric) 
        D = squareform(condensedDistances)
        
        num3Prime, num5Prime = self._geNumCellsInEachBatch()
        
        #split D
        # D is a (num3Prime + num5Prime) x (num3Prime + num5Prime)
        # the cell rows are stacked horizontally 3' on top of 5'
        # we know D[0
        D1 = D[:,0:num3Prime] # pair wise distance for ever cell and the 3 prime cells
        D2 = D[:,num5Prime:]  # pair wise distance for ever cell and the 5 prime cells
        
        our assumption that each row will contain a zero ie distance between cell_i and cell_i does not always apply
        nn1 = knn(D1) 
        nn2 = knn(D22)
        
        double check actual homework questions. should we concat or 
        keep seperate. need to computer l_k_bbknn at some pont
        bbknn = concat(nn1, nn1, axis='use rows aedwip')
        
        
        self.connectivities = None
        self.distances = None
    
    def _geNumCellsInEachBatch():
        prime3Rows = anndata.obs['Method'] == '10X_3prime'
        prime5Rows = anndata.obs['Method'] == '10X_5prime'
        num3PrimeSeries = anndata.obs.loc[prime3Rows, ['n_counts']].count()
        num5PrimeSeries = anndata.obs.loc[prime5Rows, ['n_counts']].count()
#         print("num 3':{}".format(num3PrimeSeries.values[0]))
#         print("num 5':{}".format(num5PrimeSeries.values[0]))
        return num3PrimeSeries.values[0], num5PrimeSeries.values[0]

    def get_connectivities(self,knn_indices, knn_distances):
        d, c = compute_connectivities_umap(knn_indices, 
                                          knn_distances, 
                                          knn_indices.shape[0], 
                                          knn_indices.shape[1])
        self.distances = d
        self.connectivities = c
  
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
    
# AEDWIP this is implemented in KnnG() 
#     def get_neighbors(self,D):
#
#         ''' 
#         function returns k most similar neighbors for each sample(row)
#             Input is a distance matrix calculaed from the get_distances method
#             Output is are two matrices:
#                 1. the distance of each sample (row) against every other sample (column) sorted from smallest distance (most similar) to greatest distance(least similar) for the k neighbors
#                 2. the index matrix of the sorted k nearest neighbors corresponding to the sorted distances above
#         '''
#         aedwip can we use knnG get_neighboors
#             it probably just cleaner to run pca again
#             knng = knnG(adata=None)
#             D = knng._calDistance(adataX, rep='DO NOT RUN PCA')
        
        
            
    
    def bbknn(self):
        '''
        Function computes distances for all combinations of batches.
        Fill in the below for loops
        '''
        for i in range(len(self.batch_unique)):
            #get ref batch


            for j in range(len(self.batch_unique)):

                #querry new batch with ref batch
                #this means that the distances between the querry batch and the ref_batch are calculated

                #Example:
                # if you have batches =3 and neighbors=3
                # your knn indices and distances matrices will have dims (pca_matrix.shape[0], n_batches*neighbors) = (n_obs, 9)
                # in order ot update these matrices you need to do the following:
                # for the first batch update the 0-k (n_neighbbors) columns
                # for the second batch update the k-2*k (k=n_neighbors) columns
                # for the third batch update the 2*k - 3*k (k=n_neighbors) columns
                # the indeces for the rows are always what you're querrying
                
                #the first k columns are batch1-batch1 and batch1-batch2
                #the next k columns are batch2-batch1 and batch2-batch2


        
    def l_k_bbknn(self,l=2):
        #this method makes an l subsampling of the bbknn computed in the graph() method
        # if k=4, meaning you are finding 4 neighbors between batches, and l=2:
                ## subsample 2 neighbors from batch1 and 2 neighbors from batch2
            
            
        if l >= self.neighbors_within_batch:
            raise ValueError('l cannot be equal or larger than k')
            
        self.l_knn_indices = np.zeros((self.knn_indices.shape[0], 2*l)).astype(int)
        self.l_knn_distances = np.zeros((self.knn_distances.shape[0], 2*l))
            
        
        for i in range(len(self.batch_unique)):
            #fill in loop
            pass
            for row in range(self.l_knn_indices.shape[0]):
                # select random index for each k groups
                # if k=3, self.knn_indices has 6 columns
                # so here we are going to select l random indices from first three columns, and second three columns
                # for the above example i am assuming two neighbors
                    
                #fill in loop
                pass
        pass
