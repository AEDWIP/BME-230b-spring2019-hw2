#! /usr/bin/env python


#EXAMPLE USAGE
#1st run PCA with sklearn pca module (use 50 componenets)
# euclid_knn(adata, pc=50) #make function update adata.obsm['X_pca']
# clf = knnG(adata=adata, d_metric='euclidean', n_neighbors=12) #method is always umap
# knn_indices, knn_distances = clf.get_knn()
# note that the updating of the adata object needs to be done within the get_knn() method

import logging
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import issparse
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd

# this file is strange it has functions and classes
# create a logger we can use in the functions
logger = logging.getLogger(__name__)

#pca using sklearn's pca
# this is really weird python why break out a 2 line call given we pass a string which looks like it should
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
        self.method = method
        self.reduced = None # adata.X after dimensionality has been reduced 
        
        self.get_distances(rep='pca') # this is weird python
        
        #calulcate k neighbors and umap connectivities:
        print('emptying .uns...')
        
        adata.uns['neighbors']['connectivities'] = None
        adata.uns['neighbors']['distances'] = None
        
        
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
        # this template is really weird
        
        # reduce to 50 dimensions
        tmp = None
        if rep == 'pca':
            self.reduced = myPCA(self.adata.X, 50)
            tmp = self.reduced
            
        # TODO: do not use pdist. It easier to just write a nest for loop
        # pdist might be faster if it is multi-threaded. 
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
    
    def get_neighbors(self, D):
        #fill in this method
        pass
    
    def get_umap_connectivities(self, knn_d, knn_i):
        #fill in this method
        pass
    
    def get_knn(self):
        #fill in this method
        pass
    
 
        
    