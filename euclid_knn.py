#! /usr/bin/env python


#EXAMPLE USAGE
#1st run PCA with sklearn pca module (use 50 componenets)
# euclid_knn(adata, pc=50) #make function update adata.obsm['X_pca']
# clf = knnG(adata=adata, d_metric='euclidean', n_neighbors=12) #method is always umap
# knn_indices, knn_distances = clf.get_knn()
# note that the updating of the adata object needs to be done within the get_knn() method

#from knn_to_graphModule import get_igraph_from_adjacency
import logging
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import issparse
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from scipy.interpolate.interpolate_wrapper import nearest
from sklearn.neighbors.unsupervised import NearestNeighbors

# this file is strange it has functions and classes
# create a logger we can use in the functions
logger = logging.getLogger(__name__)

# pca using sklearn's pca
# this is really weird  why break out a 2 line call given we pass a string which looks like it should
# be an enumeration that is used in a case statement? 
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
        self.distances = None # pairwise distance
        self.d_metric = d_metric
        self.n_neighbors = n_neighbors
        
        # AEDWIP: nearestNeighborsGraph is a adj matrix not a dict    
        self.nearestNeighborsGraph = dict() # graph is adj list TODO what is expected format
        self.method = method
        self.reduced = None # adata.X after dimensionality has been reduced 
        
        # wrapped initialization
        if self.adata:
            #calulcate k neighbors and umap connectivities:
            self.get_distances(rep='pca') # this is weird 
            self.get_neighbors(self.distances)
            
            print('emptying .uns...')
            
            adata.uns['neighbors']['connectivities'] = self.nearestNeighborsGraph
            adata.uns['neighbors']['distances'] = self.distances
        
        
    def update_adata(self):
        #updating adata.uns
        print('updating adata object...')
        self.adata.uns['neighbors']={}
        self.adata.uns['neighbors']['params'] = {}
        self.adata.uns['neighbors']['params']['n_neighbors']=self.n_neighbors
        self.adata.uns['neighbors']['params']['method'] = self.method
        
        self.adata.uns['neighbors']['connectivities'] = self.connectivities
        self.adata.uns['neighbors']['distances'] = self.distances
    
    def get_distances(self, rep='pca'):
        # this template is really weird. we do we need the 'rep' argument
        # are we supposed to return a distance matrix or save it to a 
        # data member
        
        # reduce to 50 dimensions
        tmp = self.adata.X
        if rep == 'pca':
            self.reduced = myPCA(self.adata.X, 50)
            tmp = self.reduced
            
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
        self.distances = squareform(condensedDistances)
        self.logger.info("self.distances.shape:{}".format(self.distances.shape))
    
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
        '''
        
        # create fast way to sort distances and look up 
        # corresponding distances
        distanceReverseIndex = { row[i]: i for i in range(len(row)) }
        
        distances = distanceReverseIndex.keys()
        # skip the first sort distance. we know it is always zero
        # it is the distance to our selves
        neighborsDistances = sorted(distances)[1: k + 1]
        
        ret = np.zeros(row.shape)
        for i in range(len(neighborsDistances)):
            distance = neighborsDistances[i]
            idx = distanceReverseIndex[ distance ]
            ret[idx] = distance
            
        return ret
    
    
    def get_neighbors(self, D):
        # aedwip how is this differ then get_knn() ???
        self.nearestNeighborsGraph = np.zeros(D.shape)
        n = D.shape[0]
        for i in range(n):
            row = D[i,:]
#             self.nearestNeighborsGraph[i] = sorted(row)[1: self.n_neighbors + 1]
            neigbors = self.findNeigborsForRow(row, self.n_neighbors)
            self.nearestNeighborsGraph[i] = neigbors
            
        return self.nearestNeighborsGraph
            
        # AEDWIP: TODO: get_igraph_from_adjacency(adjacency, directed=None):
    
    def get_umap_connectivities(self, knn_d, knn_i):
        #fill in this method
        pass
    
    def get_knn(self):
        # aedwip how is this different than get_neighbors
        #fill in this method
        pass
    
 
        
    