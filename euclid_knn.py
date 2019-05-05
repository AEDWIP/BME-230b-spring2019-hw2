#! /usr/bin/env python


#EXAMPLE USAGE
#1st run PCA with sklearn pca module (use 50 componenets)
# euclid_knn(adata, pc=50) #make function update adata.obsm['X_pca']
# clf = knnG(adata=adata, d_metric='euclidean', n_neighbors=12) #method is always umap
# knn_indices, knn_distances = clf.get_knn()
# note that the updating of the adata object needs to be done within the get_knn() method

import logging
from scipy.spatial.distance import pdist
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
        if rep == 'pca':
            self.reduced = myPCA(self.adata.X, 50)
            # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.distance.pdist.html
            self.distances = pdist(self.reduced, metric=self.d_metric)
        else :
            self.distances = pdist(self.reduced, metric=self.d_metric)
            
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
    
 
        
    