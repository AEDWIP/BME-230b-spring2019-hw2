'''
Created on May 23, 2019

@author: andrewdavidson
'''
from Cluster import Cluster
from Edge import Edge
import logging
from louvain import Louvain
from Node import Node
import numpy as np
from setupLogging import setupLogging
import unittest

class LouvainScanpyIntegrationTest(unittest.TestCase):

    setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)

    ############################################################
    def setUp(self):
        pass

    ############################################################
    def tearDown(self):
        # make sure all the logs are flushed
        # else assert we will not see partial test log entries
        logging.shutdown()

    ############################################################
    def testGetClusterAssigments(self):
        self.logger.info("BEGIN")
        
        # build lovain level tree
        listOfEdges = [(0,1), (1,0), (0,2), (2,0), (1,2), (2,1),
                        (0,3), (0,6), (0,7),
                        (3,4), (4,3), (3,5), (5,3),
                        (3,0),
                        (6,7),(7,6), (6,8), (8,6), (6,9), (9,6), (8,9), (9,8), (9,7), (7,9), 
                        (6,0), (7,0)
                       ]
        listOfWeight = [1 for i in listOfEdges]
        
        # build out init graph and calculate Q
        louvainLevel0 = Louvain.buildGraph("level0", listOfEdges, listOfWeight)
                
        # run phase I: find best cluster assignments
        louvainLevel0._phaseI(isLouvainInit=True) 
        
        ll_0_ClusterAssigments = louvainLevel0.getClusterAssigments()
        self.logger.info("level0 cluster assignments:\n{}".format(ll_0_ClusterAssigments))  
                
        # create next level and run phaseII
        # phase II consolidates clusters found in previous level
        louvainLevel1 = Louvain.buildLouvain("level1", louvainLevel0)
        louvainLevel1._calculateQ()
        
        # lets assume this is the top level.
        louvainLevel1._phaseI()     
 
        ll_1_ClusterAssigments = louvainLevel1.getClusterAssigments()
        self.logger.info("level1 cluster assignments:\n{}".format(ll_1_ClusterAssigments))  

        self.logger.info("**************** check for side effects output should be same as above")
        ll_0_ClusterAssigments = louvainLevel0.getClusterAssigments()
        self.assertEqual(ll_0_ClusterAssigments, {5: [0, 1, 2, 3, 4, 5], 9: [9, 6, 7, 8]})

        ll_1_ClusterAssigments = louvainLevel1.getClusterAssigments()
        self.assertEqual(ll_1_ClusterAssigments, {9: [9, 6, 7, 8, 0, 1, 2, 3, 4, 5]})

        self.logger.info("END\n")
           
    ############################################################
    def testRun(self):
        self.logger.info("BEGIN")
        
        # build lovain level tree
        listOfEdges = [(0,1), (1,0), (0,2), (2,0), (1,2), (2,1),
                        (0,3), (0,6), (0,7),
                        (3,4), (4,3), (3,5), (5,3),
                        (3,0),
                        (6,7),(7,6), (6,8), (8,6), (6,9), (9,6), (8,9), (9,8), (9,7), (7,9), 
                        (6,0), (7,0)
                       ]
        listOfWeight = [1 for i in listOfEdges]
        root = Louvain.run(listOfEdges, listOfWeight)
        rootClusterAssigments = root.getClusterAssigments()
        self.logger.info("louvainId: {} root cluster assigments:\n{}".format(root._louvainId, rootClusterAssigments))
        self.assertEqual(rootClusterAssigments, {9: [9, 6, 7, 8, 0, 1, 2, 3, 4, 5]})

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()