#! /usr/bin/env python

# BME-230B Spring 2019 HW 2 Question 1
# Andrew Davidson aedavids@ucsc.edu
# 
# 
# Question 1.a, 1.b see [euclid_knn.py](euclid_knn.py)
# 
# ## <span style="color:red">TODO implement 1.4</span>
# 
# ref: 
# - [ Single-Cell Analysis in Python](https://scanpy.readthedocs.io/en/stable/api/index.html#tools-tl)
# - [data exploration](exploreData.ipynb)
# - [scanpy.tl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.umap.html)
# - <span style="color:red">scanpy.api.pl no longer exists</span>
# - [scanpy.pl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.pl.umap.html#scanpy.pl.umap)

#

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

from sklearn.decomposition import PCA
import numpy as np


# create a logger for module level functions
logger = logging.getLogger(__name__)

# pca using sklearn's pca
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
    '''
    replaces scanpy pca and k nearest neighboors with our own
    
    updates:
        adata.obsm['X_pca'] 
        adata.uns['neighbors']['connectivities']
        adata.uns['neighbors']['distances']
    
    Implementation is based on template file provide with homework
    most of these functions should be private. I left the as they are becuase 
    I did not want to break the grading test harness
    
    
    usage:
        knn = KnnG(adata)
        scanpy.tl.umap(anndata)
        scanpy.pl.umap(adata)
        
    public functions
        __init__(self, adata = None, d_metric='euclidean', n_neighbors=15, method='umap')
    '''
    
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
        
        self.reduced = None # adata.X after dimensionality has been reduced 
        
        # wrapped initialization to make unit test run faster an easier easier to write
        # the functions should really be private and do not rely on data members
        # unit test should construct KnnG(adata=None)  ...
        if self.adata:
            print('emptying .uns...')
            adata.uns['neighbors']['connectivities'] = None
            adata.uns['neighbors']['distances'] = None
        
            # calulcate k neighbors and umap connectivities:
            self.D = self.get_distances(rep='pca') 
            knn_indices, knn_dist =  self.get_neighbors(self.D)
            distances,connectivities = self.get_umap_connectivities(knn_indices, knn_dist)
            self.distances = distances
            self.connectivities = connectivities

            self.update_adata()
        
    def update_adata(self):
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
        # 119,745,550 * 2 + n = 239,506,576 # coefficient of 2 give us upper and lower traingle + n is the diagonal
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
        returns a square numpy matrix. values are pair wise distances

        arguments:
            rep: representation: the algorithm to use to reduced the number of dimensions. default is 
            principle component analysis
            
            TODO: replace 'rep' with an enumeration
        '''
        self.logger.info("BEGIN")
        ret =  self._calDistance(self.adata.X, rep)
        self.logger.info("END\n")    
        return ret

    def findNeigborsForRow(self, row, k):
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
    
    def get_neighbors(self, D):
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
        knn_i = np.zeros((n,self.n_neighbors))
        knn_d = np.zeros((n,self.n_neighbors))
        for i in range(n):
            row = D[i,:]
            neigborsIdx, neigborsDist = self.findNeigborsForRow(row, self.n_neighbors)
            knn_i[i] = neigborsIdx
            knn_d[i] = neigborsDist
            
        self.logger.info("END\n")
        return knn_i, knn_d
                    
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
    
 
        
    