#! /usr/bin/env python


#EXAMPLE USAGE
#1st run PCA with sklearn pca module (use 50 componenets)
# euclid_knn(adata, pc=50) #make function update adata.obsm['X_pca']
# clf = knnG(adata=adata, d_metric='euclidean', n_neighbors=12) #method is always umap
# knn_indices, knn_distances = clf.get_knn()
# note that the updating of the adata object needs to be done within the get_knn() method

#from knn_to_graphModule import get_igraph_from_adjacency
import logging

# look into this function on scanpy's documentation to see what arguments 
# need to be passed to compute_connectivities_umap
from scanpy.neighbors import compute_connectivities_umap 

from scipy.spatial.distance import pdist, squareform
from scipy.sparse import issparse, csr_matrix

from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from scipy.interpolate.interpolate_wrapper import nearest
from sklearn.neighbors.unsupervised import NearestNeighbors
from networkx.classes.function import neighbors
from scanpy.neighbors.umap.umap_ import smooth_knn_dist


# this file is strange it has functions and classes
# create a logger we can use in the functions
logger = logging.getLogger(__name__)

# pca using sklearn's pca
# this is really weird  why break out a 3 line call given we pass a string which 
# looks like it should be an enumeration that is used in a case statement? 
def myPCA(adata, pc=15):
    '''
    input:
        adata: a numpy array
        pc: the n_components argument. should be either
            0 <= n_components <= 0
            or 
            an int greater than 0
            
            if <= 0, n_components specifies the amount of variation to preserve
            else it is the number of dimension to reduce to 
    ref: 
        - https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
    '''
    logger.info("BEGIN input shape:{} pc:{}".format(adata.shape, pc))
    
    pcaTransformer = PCA(n_components=pc)
    ret = pcaTransformer.fit_transform(adata)
        
    logger.info("END return shape:{}".format(ret.shape))
    
    return ret

    
class knnG():
    logger = logging.getLogger(__name__)

    def __init__(self, adata = None, d_metric='euclidean', n_neighbors=15, method='umap'):
        #fill in this method
        self.adata = adata 
        self.d_metric = d_metric
        self.n_neighbors = n_neighbors
        self.method = method
        
        self.D = None # the pair wise distance adjacency numpy matrix. should be very sparse
        self.connectivities = None
        self.distances = None 
        
        self.nearestNeighborsGraph = None
        self.reduced = None # adata.X after dimensionality has been reduced 
        
        # wrapped initialization to make unit test faster an easier
        # all the functionality is implemented as public functions that do
        # not rely on the data members. 
        # no need to init adata spend a lot of time recomputing distance, neighbors, ...
        if self.adata:
            print('emptying .uns...')
            adata.uns['neighbors']['connectivities'] = None
            adata.uns['neighbors']['distances'] = None
        
            # calulcate k neighbors and umap connectivities:
            self.D = self.get_distances(rep='pca') # this is weird 
            knn_indices, knn_dist =  self.get_neighbors(self.D)
#             nearestNeighborsAdjMatrix =  self.get_neighbors(self.D)
#                         
#             # convert to a sparse matrix
#             sparceNN = csr_matrix(nearestNeighborsAdjMatrix)
#             
#             # pick out the i,j tuples for non-zero value
#             rowIdx, colIdx = sparceNN.nonzero()
#             knn_i = [i for i in zip(rowIdx, colIdx)]
# 
#             # fetch the non zero distance
#             knn_d = sparceNN.data
            
            # AEDWIP: do we need to call get_umap_connectivities
            # what is going on? we already have the distances and edge tuples already 
            distances,connectivities = self.get_umap_connectivities(knn_indices, knn_dist)
            self.distances = distances
            self.connectivities = connectivities

            self.update_adata()
        
    def update_adata(self):
        #updating adata.uns
        self.logger.info('BEGIN')
        self.adata.uns['neighbors']={}
        self.adata.uns['neighbors']['params'] = {}
        self.adata.uns['neighbors']['params']['n_neighbors']=self.n_neighbors
        self.adata.uns['neighbors']['params']['method'] = self.method
        
        self.adata.uns['neighbors']['connectivities'] = self.connectivities
        self.adata.uns['neighbors']['distances'] = self.distances
        
        self.logger.info('END\n')

    
    def _calDistance(self, adataX, rep='pca'):  
        '''
        broke this out so that we can write unit tests with out having to
        run all the computation. 
        '''
        self.logger.info("BEGIN")
        # reduce to 50 dimensions
        if rep == 'pca':
            self.reduced = myPCA(adataX, 50)
            # this is really bad code
            if self.adata:
                self.adata.obsm['X_pca'] = self.reduced
            tmp = self.reduced
            self.logger.info("reduced.shape:{}".format(self.reduced.shape))
            
        #
        # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.distance.pdist.html
        # pdist() is strange
        #  X is a m x n. where m = number of examples and n is the number of dimensions
        # pdist() returns a condensed distance matrix of only the upper triangular values - the diagonal 
        # upper and lower are sysmetric, diagonal is zero
        #
        # we are told our distance matrix should be n x n where n =15,476
        # pdist() returns a shape shape: (119745550,)
        # 15,476 * 15,476  = 239,506,576
        # 119,745,550 * 2 + n = 239,506,576 # coificient of 2 give us upper and lower traingle +n is the diagonal
        #
        condensedDistances = pdist(tmp, metric=self.d_metric)
            
        # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.distance.squareform.html#scipy.spatial.distance.squareform
        # convert Converts a vector-form distance vector to a square-form distance matrix, and vice-versa
        ret = squareform(condensedDistances)
        self.logger.info("self.distances.shape:{}".format(ret.shape))
        
        self.logger.info("END\n")
        return ret
    
    def get_distances(self, rep='pca'):
        '''
        AEDWIP: TODO: #returns a square numpy "adjacency" matrix. values are pair wise distances
        '''
        # this template is really weird. we do we need the 'rep' argument
        # are we supposed to return a distance matrix or save it to a 
        # data member
        self.logger.info("BEGIN")
        ret =  self._calDistance(self.adata.X, rep)
        self.logger.info("END\n")    
        return ret
    
#     def get_neighbors(self, D):
#         # aedwip how is this differ then get_knn() ???
#         n = D.shape[0]
#         for i in range(n):
#             row = D[i,:]
#             self.nearestNeighborsGraph[i] = sorted(row)[1: self.n_neighbors + 1]
#             
#         # AEDWIP: TODO: get_igraph_from_adjacency(adjacency, directed=None):


    def findNeigborsForRow(self, row, k):
        '''
        arguments:
            row: a row in the pair wise distance matrix
            k: the number of nearest neighbors to find
        
        returns:
            a row in adjacency matrix format of k nearest neighbors
            two numpy arrays. The values in the first array are the indices 
            for the row argument nearest neighbors. The values of the second
            array are the distances
        '''
        
        # create fast way to sort distances and look up 
        # corresponding distances
        distanceReverseIndex = { row[i]: i for i in range(len(row)) }
        
        distances = distanceReverseIndex.keys()
        # skip the first sort distance. we know it is always zero
        # it is the distance to our selves
        neighborsDistances = sorted(distances)[1: k + 1]
        
        retIdx = np.zeros(k)
        retDist = np.zeros(k)
        for i in range(len(neighborsDistances)):
            distance = neighborsDistances[i]
            idx = distanceReverseIndex[ distance ]
            retIdx[i] = idx
            retDist[i] = distance
            
        return retIdx,retDist
    
    def get_neighbors(self, D):
        '''
        arguments
            D: pairwise distance matrix
            
        returns and adjacency numpy matrix. It should be very sparse
        
        AEDWIP: TODO: should we return a sparse matrix
        '''
        self.logger.info("BEGIN")
        #nearestNeighborsAdjMatrix = np.zeros(D.shape)
        n = D.shape[0]
        knn_i = np.zeros((n,self.n_neighbors))
        knn_d = np.zeros((n,self.n_neighbors))
        for i in range(n):
            row = D[i,:]
#             neigbors = self.findNeigborsForRow(row, self.n_neighbors)
#             nearestNeighborsAdjMatrix[i] = neigbors
            neigborsIdx, neigborsDist = self.findNeigborsForRow(row, self.n_neighbors)
            knn_i[i] = neigborsIdx
            knn_d[i] = neigborsDist
            
        self.logger.info("END\n")
        #return nearestNeighborsAdjMatrix
        return knn_i, knn_d
            
        # AEDWIP: TODO: get_igraph_from_adjacency(adjacency, directed=None):
        # looks like the adjMatrix should be in sparce format
        
    def get_umap_connectivities(self, knn_indices, knn_dists ):
        '''
        ref:
            https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.Neighbors.html?highlight=scanpy.neighbors
        '''
        self.logger.info("BEGIN")
        # you have to read the code to figure out how to call compute_connectivities_umap
        n_obs = self.adata.X.shape[0]
        
        self.logger.info("type(knn_dists):{} knn_dists.shape:{}".format(type(knn_dists), knn_dists.shape))
        self.logger.info("knn_dists[0:3:\n{}".format(knn_dists[0:3]))

        self.logger.info("type(knn_indices):{}".format(type(knn_indices)))
        # does not have 'shape' self.logger.info("knn_indices.shape:{}".format(knn_indices.shape)
        self.logger.info("knn_indices[0:3]:\n{}".format(knn_indices[0:3]))
        # ??? 
        
        try:
#             knn_indices should be our nearestNeighboes adj matrx
#             knn_dists is same struct
            
#             knn_indices is ?? array or rows. each contains the idx of its nearest neighbors
#             knn_dists is a similar struct, the values are the distance to the neighboors  
            distances,connectivities =  compute_connectivities_umap(knn_indices, 
                                                                knn_dists,
                                                                n_obs, 
                                                                self.n_neighbors)
        except Exception as e: # catch *all* exceptions
            # problem logger does not always flush
            self.logger.info(e)
            self.logger.error(e)
            raise e
        
        self.logger.info("type(distances):{}".format(type(distances)))
        self.logger.info("distances:\n{}".format(distances))
        
        self.logger.info("type(connectivities):{}".format(type(connectivities)))
        self.logger.info("connectivities:\n{}".format(connectivities))

        self.logger.info("END")
        
        return distances,connectivities
    
#     def get_knn(self):
#         '''
#         returns:
#             knn_indices 
#             knn_distances 
#         '''
#         self.logger.info("BEGIN")
#         # AEWIP where / when is this function called?
#         
#         # this is weird everything should be initialiized in __init__
#         self.update_adata()
#         self.logger.info("END")
# 
#         return self.connectivities, self.distances
    
 
        
    