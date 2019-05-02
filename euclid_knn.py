#! /usr/bin/env python


#EXAMPLE USAGE
#1st run PCA with sklearn pca module (use 50 componenets)
# euclid_knn(adata, pc=50) #make function update adata.obsm['X_pca']
# clf = knnG(adata=adata, d_metric='euclidean', n_neighbors=12) #method is always umap
# knn_indices, knn_distances = clf.get_knn()
# note that the updating of the adata object needs to be done within the get_knn() method



from scipy.sparse import issparse
import numpy as np
import pandas as pd

#pca using sklearn's pca
def pca(adata, pc=15):
    pass
    

    
class knnG():
    def __init__(self, adata = None, d_metric='euclidean', n_neighbors=15, method='umap'):
        #fill in this method
        
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
        #fill in this method
        pass
    
    def get_neighbors(self, D):
        #fill in this method
        pass
    
    def get_umap_connectivities(self, knn_d, knn_i):
        #fill in this method
        pass
    
    def get_knn(self):
        #fill in this method
        pass
    
 
        
    