'''
Created on May 5, 2019

@author: andrewdavidson
'''

from euclid_knn import knnG
#from knn_to_graphModule import get_igraph_from_adjacency
import logging
import numpy as np
import scanpy.api as sc
from scipy.spatial.distance import pdist
from setupLogging import setupLogging
import unittest

class euclid_knnTest(unittest.TestCase):
    setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)

    def setUp(self):
        pass


    def tearDown(self):
        # make sure all the logs are flushed
        # if assert we will not see partial test log entries
        logging.shutdown()

    ############################################################
    def testPairWise(self):
        '''
        pdist assumes row major
        X is m x n. 
        m examples in n dimensions
        
        pdist() return is a condensed vector. see comments in euclid_knn 
        get_distance() for the details
        '''
        self.logger.info("BEGIN ")
        
        X = np.array([[1, 2], 
                      [3, 4],
                      [5, 6]])
                
        expected = np.array([
            np.sqrt( (1-3)**2 + (2-4)**2 ),
            np.sqrt( (1-5)**2 + (2-6)**2 ),
            np.sqrt( (3-5)**2 + (4-6)**2 )
            ])
        
        ret = pdist(X, metric='euclidean')
        self.logger.info("x:\n{}".format(X))
        self.logger.info("ret:\n{}".format(ret))
        self.logger.info("expected:\n{}".format(expected))
        
        self.assertTrue( (ret == expected).all() )
        
        self.logger.info("END \n")

    ############################################################
    def testGetDistance(self):
        '''
        TODO:
        '''
        self.logger.info("BEGIN ")
        
        anndata = sc.read("../PBMC.merged.h5ad")
        knng = knnG(anndata)
        
        self.logger.info("knng.reduced.shape: {}".format(knng.reduced.shape))
        self.logger.info("knng.distances.shape: {}".format(knng.distances.shape))
        
        self.assertTrue( (knng.reduced.shape   == (15476, 50)) )
        self.assertTrue( (knng.distances.shape == (15476, 15476)) )
        
        print("knng.distna:\n{}".format(knng.distances[0:5,0:5]))
        
        self.logger.info("END \n")

    ############################################################
    @unittest.skip("skip this test it assume adj list not adj matrix")
    def testGetNeighborsAdjList(self):
        '''
        TODO:
        '''
        self.logger.info("BEGIN ")
        
        anndata = sc.read("../PBMC.merged.h5ad")
        knng = knnG(anndata)
        
        # call get_neighbors using our test data
        D = np.array([
                        [ 0, 1, 2, 3, 4, 5],
                        [11, 0, 9, 8, 7, 6],
                        [12, 3, 0, 5, 4, 2]
                    ])
        
        knng.n_neighbors = 4
        knng.nearestNeighborsGraph = {}
        knng.get_neighbors(D)
        self.logger.info("nearestNeighbors:\n{}".format(knng.nearestNeighborsGraph))
        
        expected = {0: [1, 2, 3, 4], 1: [6, 7, 8, 9], 2: [2, 3, 4, 5]}
        self.assertTrue(knng.nearestNeighborsGraph == expected)
        
        self.logger.info("END \n")

      
    ############################################################
    def testGetNeighborsAdjMatrix(self):
        '''
        TODO:
        '''
        self.logger.info("BEGIN ")
        
        mockDistance = np.array([
                        [ 0, .1, .2, .3, .4, .5],
                        [.11, 0, .9, .8, .7, .6],
                        [.12, .3, 0, .5, .4, .2]
                    ])
        
        knng = knnG(None)
        
        knng.n_neighbors = 4
        knng.nearestNeighborsGraph = {}
        ret = knng.get_neighbors(mockDistance)
        self.logger.info("nearestNeighbors:\n{}".format(ret))

        expected = np.array([
                            [0.,   0.1,  0.2,  0.3,  0.4,  0.  ],
                             [0.11, 0.,   0.,   0.8,  0.7,  0.6 ],
                             [0.12, 0.3,  0.,   0.,   0.4,  0.2 ]
                             ])
        
        self.logger.info("expected:\n{}".format(expected))
        
        self.assertTrue( (ret == expected).all() )
        
        self.logger.info("END \n")
        
    ############################################################        
    def test_get_igraph_from_adjacency(self):
        self.logger.info("BEGIN ")
        
        adjacency = {
                1:[1, 2, 3],
                2:[2, 3],
                3:[4]
            }
        ret = get_igraph_from_adjacency(adjacency)
        self.logger.info("ret:\n{}".format(ret))
        self.logger.info("END \n ")



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()