'''
Created on May 25, 2019

@author: andrewdavidson
'''
import logging
from louvain import Louvain
import numpy as np
# import scanpy.api as sc
# from scipy.spatial.distance import pdist
from setupLogging import setupLogging
import unittest
from Edge import Edge
from Node import Node
from Cluster import Cluster

############################################################
class testDisjointGraphs(unittest.TestCase):

    setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)

    ############################################################
    def setUp(self):
        pass

    ############################################################
    def tearDown(self):
        # make sure all the logs are flushed
        # else if assert we will not see partial test log entries
        logging.shutdown()
        
    ############################################################
    def testDJGraph(self):
        self.logger.info("BEGIN")
        
        listOfEdges = [(0,1), (1,0), (2,3), (3,2) ]
        
        # all other test assume weights are 1. set to a big number
        # that will make it easy to spot bugs in summary stats.
        listOfWeight = [5, 5, 10, 10]   
        numRows = 4 # number of nodes
        
        louvainLevel0 = Louvain.buildGraph("level 0", listOfEdges, listOfWeight)
        # check Q
        louvainLevel0._calculateQ()
        self.logger.info("after  buildGraph() louvainLevel0._Q:{}".format(louvainLevel0._Q))
        self.assertAlmostEqual(louvainLevel0._Q, 0)        
        
        # check is initialization of graph correct
        print()
        for cluster in louvainLevel0._clustersLookup.values():
            self.logger.info("{}".format(cluster))
           
        print() 
        for node in louvainLevel0._nodeLookup.values():
            self.logger.info("\nnode {}".format(node))
            print()
            
        expectedNodeConfigDict = {
            0:"clusterId:0 nodeId:0 numEdges:1 adjEdgeWeights:10 _weightsInClusterDict{1: 10} _edgesDict{1: srcId:0 targetId:1 weight:10}",
            1: "clusterId:1 nodeId:1 numEdges:1 adjEdgeWeights:10 _weightsInClusterDict{0: 10} _edgesDict{0: srcId:1 targetId:0 weight:10}",
            2:"clusterId:2 nodeId:2 numEdges:1 adjEdgeWeights:10 _weightsInClusterDict{3: 10} _edgesDict{3: srcId:2 targetId:3 weight:10}",
            3:"nodeId:3 numEdges:1 adjEdgeWeights:10 _weightsInClusterDict{2: 10} _edgesDict{2: srcId:3 targetId:2 weight:10}"
        }
            
        # check kiin
        for node in louvainLevel0._nodeLookup.values():
            self.logger.info("nodeId:{} _weightsInClusterDict:{}".format(node._nodeId, node._weightsInClusterDict))
            
       
        # run phase I
        louvainLevel0._phaseI(numRows, isLouvainInit=True) # TODO: can probably get rid of isLouvainInit
        self.logger.info("after phase I() louvainLevel0:\n{}".format(louvainLevel0))
        l0Assignments = louvainLevel0.getClusterAssigments()
        self.logger.info("l0Assigments cluster assignments:\n{}".format(l0Assignments))        
        
        # check Q
        louvainLevel0._calculateQ()
        self.logger.info("after phase I   louvainLevel0._Q:{}".format(louvainLevel0._Q))
        self.assertAlmostEqual(louvainLevel0._Q, 0.7222222222222221)
             
        # build next level
        louvainLevel1 = Louvain.buildLouvain("level 1 ", louvainLevel0)
        self.logger.info("after buildLouvain louvainLevel1\n{}".format(louvainLevel1))
        
        # phase II
        louvainLevel1._phaseII(isLouvainInit=False) # TODO: can probably get rid of isLouvainInit)
        self.logger.info("after phaseII() louvainLevel1  this log line looks funnny:\n{}".format(louvainLevel1)) 
        l1Assignments = louvainLevel1.getClusterAssigments()
        self.logger.info("louvainLevel1 cluster assignments:\n{}".format(l1Assignments))
     
        self.logger.info("END\n")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()