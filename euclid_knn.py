#! /usr/bin/env python

# BME-230B Spring 2019 HW 2 Question 1
# James Casaletto, Andrew Davidson, Yuanqing Xue, Jim Zheng
# 
# 
# Question 1.a, 1.b see [euclid_knn.py](euclid_knn.py)
# 
# 
# ref: 
# - [ Single-Cell Analysis in Python](https://scanpy.readthedocs.io/en/stable/api/index.html#tools-tl)
# - [data exploration](exploreData.ipynb)
# - [scanpy.tl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.umap.html)
# - <span style="color:red">scanpy.api.pl no longer exists</span>
# - [scanpy.pl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.pl.umap.html#scanpy.pl.umap)

import logging
from scanpy.neighbors import compute_connectivities_umap 
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
import numpy as np

################################################################################
class KnnG():
    '''
    replaces scanpy pca and k nearest  neighbours  with our own
    
    updates:
        adata.obsm['X_pca'] 
        adata.uns['neighbors']['connectivities']
        adata.uns['neighbors']['distances']
    
    usage:
        knn = KnnG(adata)
        scanpy.tl.umap(anndata)
        scanpy.pl.umap(adata)
        
    public functions
        __init__(self, adata = None, d_metric='euclidean', n_neighbors=15, method='umap',
                    runPCA=True, nPC=50)
    '''
    
    logger = logging.getLogger(__name__)

    ######################################################################
    def __init__(self, adata = None, d_metric='euclidean', n_neighbors=15, 
                 method='umap', runPCA=True, nPC=50):
        '''
        arguments:
            adata: unit test should pass None and init members as needed
            nPC: if runPCA == True should be the number of principle components you want
                to reduce adata.X to 
        '''
        self._adata = adata 
        self._d_metric = d_metric
        self._n_neighbors = n_neighbors
        self._method = method
        
        self._D = None # the pair wise distance a numpy matrix. 
        self._connectivities = None
        self._distances = None 
        
        # wrapped initialization to make unit test run faster an easier to write
        # unit test should construct KnnG(adata=None)
        if self._adata:
            print('emptying .uns...')
            self._adata.uns['neighbors']['connectivities'] = None
            self._adata.uns['neighbors']['distances'] = None
        
            # calculate k neighbors and umap connectivities:
            if runPCA:
                self._PCA(nPC)
            
            self._calDistance()
            knn_indices, knn_dist =  self._get_neighbors(self.D)
            self._get_umap_connectivities(knn_indices, knn_dist)

            self._update_adata()
        

    ######################################################################    
    def _calDistance(self):  
        '''
        calculates pair wise distances. results store in self.D
        
        input:
        '''
        self.logger.info("BEGIN")
            
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
        # 119,745,550 * 2 + n = 239,506,576 # coefficient of 2 give us upper and lower traingle + n is the diagonal
        #
        condensedDistances = pdist(self._adata.obsm['X_pca'], metric=self._d_metric)
            
        # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.distance.squareform.html#scipy.spatial.distance.squareform
        # convert Converts a vector-form distance vector to a square-form distance matrix, and vice-versa
        self._D = squareform(condensedDistances)
        self.logger.info("self.distances.shape:{}".format(self._D.shape))
        
        self.logger.info("END\n")
    
    ######################################################################    
    def _findNeigborsForRow(self, row, k):
        '''
        arguments:
            row: a row in the pair wise distance matrix
            k: the number of nearest neighbors to find
        
        returns:
            a row in adjacency matrix format of k nearest neighbors as
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
    
    ######################################################################
    def _get_neighbors(self, D):
        '''
        arguments
            D: pairwise distance matrix
            
        returns:
            a row in adjacency matrix format of k nearest neighbors as
            two numpy arrays. The values in the first array are the indices 
            for the row argument nearest neighbors. The values of the second
            array are the distances        
        '''
        self.logger.info("BEGIN")
        n = D.shape[0]
        knn_i = np.zeros((n,self._n_neighbors))
        knn_d = np.zeros((n,self._n_neighbors))
        for i in range(n):
            row = D[i,:]
            neigborsIdx, neigborsDist = self._findNeigborsForRow(row, self._n_neighbors)
            knn_i[i] = neigborsIdx
            knn_d[i] = neigborsDist
            
        self.logger.info("END\n")
        return knn_i, knn_d
      
    ######################################################################              
    def _get_umap_connectivities(self, knn_indices, knn_dists ):
        '''
        ref:
            https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.Neighbors.html?highlight=scanpy.neighbors
            
        results are store in
            self.distances = distances
            self.connectivities = connectivities
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
        
        # knn_indices is ?? array of rows. each contains the idx of its nearest neighbors
        # knn_dists is a similar structure, the values are the distance to the neighboors  
        distances,connectivities =  compute_connectivities_umap(knn_indices, 
                                                            knn_dists,
                                                            n_obs, 
                                                            self.n_neighbors)
        
        self.logger.info("type(distances):{}".format(type(distances)))
        self.logger.info("distances:\n{}".format(distances))
        
        self.logger.info("type(connectivities):{}".format(type(connectivities)))
        self.logger.info("connectivities:\n{}".format(connectivities))

        self.logger.info("END")
        
        self.distances = distances
        self.connectivities = connectivities
    
    ######################################################################
    def _PCA(self, npc=15):
        '''
        input:
            npc: the n_components argument. should be either
                0 <= n_components <= 0
                or 
                an int greater than 0
                
                if <= 0, n_components specifies the amount of variation to preserve
                else it is the number of dimension to reduce to 
                
        results are stored in self.adata.obsm['X_pca']
        ref: 
            - https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
        '''
        self.logger.info("BEGIN input shape:{} pc:{}".format(self._adata.shape, npc))
        
        pcaTransformer = PCA(n_components=npc)
        self._adata.obsm['X_pca'] = pcaTransformer.fit_transform(self._adata.X)

        self.logger.info("END return shape:{}".format(self._adata.obsm['X_pca'].shape)) 

    ######################################################################
    def _update_adata(self):
        '''
        The trick to replacing scanpy implementation with our own is to 
        update the anndata object with our intermediate values
        
        we replace the scanpy version of PCA by updating
        adata.obsm['X_pca'] = our PCA() output
        
        we replace the scanpy version of k-nearest-neighbors by updating
        self.adata.uns['neighbors']['connectivities'] = our knn() output
        self.adata.uns['neighbors']['distances'] = out knn() output
        '''
        self.logger.info('BEGIN')
        self.adata.uns['neighbors']={}
        self.adata.uns['neighbors']['params'] = {}
        self.adata.uns['neighbors']['params']['n_neighbors']=self.n_neighbors
        self.adata.uns['neighbors']['params']['method'] = self.method
        
        self.adata.uns['neighbors']['connectivities'] = self.connectivities
        self.adata.uns['neighbors']['distances'] = self.distances
        
        self.logger.info('END\n')
