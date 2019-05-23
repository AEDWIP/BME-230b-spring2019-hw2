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
    def checkClusters(self, expected, clustersDict):
        self.logger.info("BEGIN")
        for clusterId, cluster in clustersDict.items():
            self.logger.info("clusterId:{}".format(clusterId))
            self.assertEqual(len(cluster._nodeList), expected[clusterId]['numNodes'])
            self.assertEqual(cluster._weightsInsideCluster, expected[clusterId]['weightsInsideCluster'])
            self.assertEqual(cluster._totalWeight, expected[clusterId]['totalWeight']) 
            print()
                   
        self.logger.info("END\n")
        
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
            
        self.checkClusters(expected, louvainLevel0._clusters)
            
        self.assertEqual(louvainLevel0._Q, 0.5)
        
        self.logger.info("END\n")
        
    ############################################################
    def testPhaseI(self):
        self.logger.info("BEGIN")
        
        listOfEdges = [(0,1), (1,0),    (0,2), (2,0),    (0,3), (3,0),
                       (1,2), (2,1),    (3,4), (4,3) ]
        listOfWeight = [1 for i in listOfEdges]
        louvainLevel0 = Louvain.buildGraph("testPhaseI graph", listOfEdges, listOfWeight)
        
        louvainLevel0._phaseI(isLouvainInit=True) 
        
        for clusterId, cluster in louvainLevel0._clusters.items():
            print('')
            self.logger.info("cluserId:{}".format(clusterId))
            self.logger.info(cluster)
            
        expected = {
            0 : {'clusterId':0, 'numNodes':0, 'weightsInsideCluster':0, 'totalWeight':0},
            1 : {'clusterId':1, 'numNodes':0, 'weightsInsideCluster':0, 'totalWeight':0},
            2 : {'clusterId':2, 'numNodes':0, 'weightsInsideCluster':0, 'totalWeight':0},
            3 : {'clusterId':3, 'numNodes':3, 'weightsInsideCluster':8, 'totalWeight':7},
            4 : {'clusterId':4, 'numNodes':2, 'weightsInsideCluster':2, 'totalWeight':3}
            }
        
        self.checkClusters(expected, louvainLevel0._clusters)

        c3 = louvainLevel0._clusters[3]
        for n in c3._nodeList:
            self.logger.info("clusterId:{}, nodeId:{}".format(c3._clusterId, n._nodeId))
            
        print('')
        c4 = louvainLevel0._clusters[4]
        for n in c4._nodeList:
            self.logger.info("clusterId:{}, nodeId:{}".format(c4._clusterId, n._nodeId))  
            
        expectedNodesInCluster = {
            3: [0,1,2],
            4: [3,4]
            }
        
        for clusterId in expectedNodesInCluster.keys():
            c = louvainLevel0._clusters[clusterId]  
            nodeIds = [n._nodeId for n in c._nodeList]     
            self.assertEqual(sorted(nodeIds), sorted(expectedNodesInCluster[clusterId]))
        
        self.logger.info("END\n")        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()