'''
Created on May 14, 2019

@author: andrewdavidson
'''
import logging
import numpy as np
from setupLogging import setupLogging

import unittest
from euclid_bbknn import bbknn_graph



class Test(unittest.TestCase):
    setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)

    def setUp(self):
        # make test reproducable
        np.random.seed(42)

    def tearDown(self):
        # make sure all the logs are flushed
        # if assert we will not see partial test log entries
        logging.shutdown()

    def testl_k_bbknn(self):
        self.logger.info("BEGIN")

        # two batches
        # batch0 has 3 cells
        # batch1 has 3 cells
        pairwiseDist=np. array([
                                [2,3,4,6,5,4,3],
                                [12,11,13,16,15,14,13],
                                [22,21,23,26,25,24,23],
                                [32,31,33,36,35,34,33],
                                [42,41,43,46,45,44,43],
                                [52,51,53,56,55,54,53],
                                [62,61,63,66,65,64,63]
                                ])
        
        bb2nnIdx = np.array([
                                [1, 2, 6, 4],
                                [1, 2, 6, 4],
                                [1, 2, 6, 4],
                                [1, 2, 6, 4],
                                [1, 2, 6, 4],
                                [1, 2, 6, 4],
                                [1, 2, 6, 4]
            ])
        
        bb2nnDists = np.array([
                                    [ 1,  2,  3,  4],
                                    [11, 12, 13, 14],                                    
                                    [21, 22, 23, 24],
                                    [31, 32, 33, 34],
                                    [41, 52, 53, 54],
                                    [51, 62, 63, 64],
                                    [61, 62, 63, 64]
                                    ])
        
        bbknn = bbknn_graph(adata=None, batchlabel = None, 
                            neighbors_within_batch=2, pcs=None, method=None,
                            batch_unique=2)
        
        
        bbknn.knn_indices = bb2nnIdx
        bbknn.knn_distances = bb2nnDists
        
        bbknn.l_k_bbknn(l=1)
        
        # get results
        retl_knn_indices = bbknn.l_knn_indices
        self.logger.info("retl_knn_indices:\n{}".format(retl_knn_indices))
        
        retl_knn_distances = bbknn.l_knn_distances
        self.logger.info("retl_knn_distances:\n{}".format(retl_knn_distances))

        expectedIdx = np.array([[1, 4],[1, 6],[1, 4],
                                [1, 6],[1, 4],[1, 6],[1, 6]])
        np.testing.assert_array_equal(expectedIdx, retl_knn_indices)
        
        expectedDist = np.array([[ 1.,4.],[11., 13.],[21., 24.],
                                 [31., 33.],[41., 54.],[51., 63.]])
        np.testing.assert_array_equal(expectedDist, retl_knn_distances)

       

        self.logger.info("END\n")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()