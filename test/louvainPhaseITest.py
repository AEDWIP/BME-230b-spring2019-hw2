'''
Created on May 22, 2019

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

############################################################
class LouvianPhaseITest(unittest.TestCase):
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
    def testSimplePhaseI(self):
        self.logger.info("BEGIN")
        
        # create a trival graph.
        listOfEdges = [(0,1), (1,2)]
        listOfWeight = [1 for i in listOfEdges]
        louvainLevel0 = Louvain.buildGraph("trival graph", listOfEdges, listOfWeight)
        
        louvainLevel0._phaseI(isLouvainInit=True) 
        
        for clusterId, cluster in louvainLevel0._clusters.items():
            print()
            self.logger.info("cluserId:{}".format(clusterId))
            self.logger.info(cluster)
            
        expected = {
                    0 : {'clusterId':0, 'numNodes':0, 'weightsInsideCluster':0, 'totalWeight':0},
                    1 : {'clusterId':1, 'numNodes':3, 'weightsInsideCluster':4, 'totalWeight':4},
                    2 : {'clusterId':2, 'numNodes':0, 'weightsInsideCluster':0, 'totalWeight':0}
                    }
        for clusterId, cluster in louvainLevel0._clusters.items():
            self.assertEqual(len(cluster._nodeList), expected[clusterId]['numNodes'])
            self.assertEqual(cluster._weightsInsideCluster, expected[clusterId]['weightsInsideCluster'])
            self.assertEqual(cluster._totalWeight, expected[clusterId]['totalWeight'])
            
        self.assertEqual(louvainLevel0._Q, 0.5)
        
        self.logger.info("END\n")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()