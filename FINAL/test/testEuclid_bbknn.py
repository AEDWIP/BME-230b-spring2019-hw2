'''
Created on May 14, 2019

@author: andrewdavidson
'''
import logging
import numpy as np
from setupLogging import setupLogging
import scanpy as sc

import unittest
from euclid_bbknn import bbknn_graph


######################################################################                
class Eculid_bbknnTest(unittest.TestCase):
    setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)
    
    ######################################################################                
    def setUp(self):
        # make test reproducable
        np.random.seed(42)

    ######################################################################                
    def tearDown(self):
        # make sure all the logs are flushed
        # if assert we will not see partial test log entries
        logging.shutdown()

    ######################################################################                
    def test_l_k_bbknn(self):
        self.logger.info("BEGIN")

        # two batches
        # batch0 has 3 cells
        # batch1 has 3 cells
        
#         pairwiseDist=np. array([
#                                 [2,3,4,6,5,4,3],
#                                 [12,11,13,16,15,14,13],
#                                 [22,21,23,26,25,24,23],
#                                 [32,31,33,36,35,34,33],
#                                 [42,41,43,46,45,44,43],
#                                 [52,51,53,56,55,54,53],
#                                 [62,61,63,66,65,64,63]
#                                 ])
        
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
                                    [41, 42, 43, 44],
                                    [51, 52, 53, 54],
                                    [61, 62, 63, 64]
                                    ])
        
        bbknn = bbknn_graph(adata=None, 
                            neighbors_within_batch=2, pcs=None, method=None)
        
        # init bbknn with our test data
        bbknn._knn_indices = bb2nnIdx
        bbknn._knn_distances = bb2nnDists
        bbknn._numBatches = 2
        
        bbknn._l_k_bbknnImplementation(l=1)
        
        # get results
        retl_knn_indices = bbknn._l_knn_indices
        self.logger.info("retl_knn_indices:\n{}".format(retl_knn_indices))
        
        retl_knn_distances = bbknn._l_knn_distances
        self.logger.info("retl_knn_distances:\n{}".format(retl_knn_distances))

        expectedIdx = np.array([[1, 4],
                                [1, 6],
                                [1, 4],
                                [1, 6],
                                [1, 4],
                                [1, 6],
                                [1, 6]])
        np.testing.assert_array_equal(expectedIdx, retl_knn_indices)
        
        expectedDist = np.array([[ 1.,4.],
                                 [11., 13.],
                                 [21., 24.],
                                 [31., 33.],
                                 [41., 44.],
                                 [51., 53.], 
                                 [61., 63.]])
        np.testing.assert_array_equal(expectedDist, retl_knn_distances)
        
        # make sure random sub setting works as expected
        bbknn._l_k_bbknnImplementation(l=1)
        retl_knn_indices2 = bbknn._l_knn_indices
        self.logger.info("retl_knn_indices2:\n{}".format(retl_knn_indices2))
        
        expectedIdx2 = np.array([[2, 6],
                                 [2, 4],
                                 [2, 6],
                                 [2, 6],
                                 [2, 4],
                                 [2, 4],
                                 [2, 4]])
        
        np.testing.assert_array_equal(expectedIdx2, retl_knn_indices2)
        
        
        retl_knn_distances2 = bbknn._l_knn_distances
        self.logger.info("retl_knn_distances2:\n{}".format(retl_knn_distances2))
        expectedDist2 = np.array([[ 2.,  3.],
                                  [12., 14.],
                                  [22., 23.],
                                  [32., 33.],
                                  [42., 44.],
                                  [52., 54.],
                                  [62., 64.]])
        
        np.testing.assert_array_equal(expectedDist2, retl_knn_distances2)


        self.logger.info("END\n")
       
   ######################################################################                 
    def testBbknn(self):
        self.logger.info("BEGIN")
        # two batches
        # batch0 has 3 cells
        # batch1 has 4 cells
        pairwiseDist=np. array([
                                [ 0,  1,  3,  6,  5,  4,  3],
                                [12,  0, 13, 16, 15, 14, 13],
                                [22, 21,  0, 26, 25, 24, 23],
                                [32, 31, 33,  0, 35, 34, 33],
                                [42, 41, 43, 46,  0, 44, 43],
                                [52, 51, 53, 56, 55,  0, 53],
                                [62, 61, 63, 66, 65, 64,  0]
                                ])
        
        expectedBB2NNIdx = np.array([
                            [1, 2, 6, 5],
                            [0, 2, 6, 5],
                            [1, 0, 6, 5],
                            [1, 0, 6, 5],
                            [1, 0, 6, 5],
                            [1, 0, 6, 4],
                            [1, 0, 5, 4]
                            ], dtype=float)
        
        expectedBB2NNDists = np.array([
                                [ 1,  3,  3,  4],
                                [12, 13, 13, 14],                                    
                                [21, 22, 23, 24],
                                [31, 32, 33, 34],
                                [41, 42, 43, 44],
                                [51, 52, 53, 55],
                                [61, 62, 64, 65]
                                ])      
        
        
        bbknn = bbknn_graph(None, neighbors_within_batch=2)  
        
        batchCounts= [('0', 3), ('1', 4)]
        retBBKNNIdx,retBBKNNDist = bbknn._bbknn(D=pairwiseDist, batchCounts=batchCounts)
        self.logger.info("retBBKNNIdx:\n{}".format(retBBKNNIdx))        
        self.logger.info("expectedBB2NNIdx:\n{}".format(expectedBB2NNIdx))
        np.testing.assert_array_equal(expectedBB2NNIdx, retBBKNNIdx)

        self.logger.info("retBBKNNDist:\n{}".format(retBBKNNDist))
        self.logger.info("expectedBB2NNDists:\n{}".format(expectedBB2NNDists))
        np.testing.assert_array_equal(expectedBB2NNDists, retBBKNNDist)


        self.logger.info("END\n")
        
   ######################################################################                 
    def testCalcSplits(self):
        self.logger.info("BEGIN")
        
        bbknn = bbknn_graph(None)
        batchCounts= [('0', 3), ('1', 4)]
        splitsLocations = bbknn._calcSplits(batchCounts)
        self.assertEqual(splitsLocations, [3,7])
        self.logger.info("splitsLocations:{}".format(splitsLocations))
        
        D= np. array([
                    [ 0,  1,  3,  6,  5,  4,  3],
                    [12,  0, 13, 16, 15, 14, 13],
                    [22, 21,  0, 26, 25, 24, 23],
                    [32, 31, 33,  0, 35, 34, 33],
                    [42, 41, 43, 46,  0, 44, 43],
                    [52, 51, 53, 56, 55,  0, 53],
                    [62, 61, 63, 66, 65, 64,  0]
                    ])
        self.logger.info("D.shape:{}".format(D.shape))
                
        byCols = 1
        splits = np.split(D, splitsLocations, axis=byCols)
#         splits = np.split(D, [3,:], axis=byCols)
        self.logger.info("AEDWIP len(splits):{}".format(len(splits)))
        
        # np.split(D, [3,6] returns
        # [:3], [3:6], [6:]
        # we need to remove this last split
        del splits[-1]
        for split in splits:
            self.logger.info("AEDWIP split.shape():{}".format(split.shape))
            self.logger.info("split\n{}\n".format(split))
            
        self.assertEqual(len(splits), 2)
        
        self.logger.info("END\n")

   ######################################################################
    @unittest.skip("skip this test it takes over 3 mins")                 
    def test_l_k_adata(self):
        self.logger.info("BEGIN")
        
        anndata = sc.read("../PBMC.merged.h5ad")
        bb6nn = bbknn_graph(anndata, neighbors_within_batch=6, runPCA=False, pcs=50)
        bb6nn.l_k_bbknn(l=3)
        
        l_k_d  = anndata.uns['neighbors']['distances']
        self.logger.info("type(l_k_d):{}".format(type(l_k_d)))
        self.logger.info("l_k_d.shape:{}".format(l_k_d.shape))
        self.logger.info("l_k_d:{}\n".format(l_k_d))

        self.logger.info("END\n")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()