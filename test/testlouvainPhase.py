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
class LouvianPhaseTest(unittest.TestCase):
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
            msg="clusterId:{}".format(clusterId)
            self.assertEqual(len(cluster._nodeList), expected[clusterId]['numNodes'], msg)
            self.assertEqual(cluster._weightsInsideCluster, expected[clusterId]['weightsInsideCluster'], msg)
            self.assertEqual(cluster._totalWeight, expected[clusterId]['totalWeight'], msg) 
                   
        self.logger.info("END\n")
        
    ############################################################
    def testSimplePhaseI(self):
        self.logger.info("BEGIN")
        
        # create a trival graph.
        listOfEdges = [(0,1), (1,2)]
        listOfWeight = [1 for i in listOfEdges]
        louvainLevel0 = Louvain.buildGraph("trival graph", listOfEdges, listOfWeight)
        louvainLevel0._calculateQ()
        numRows = 3 # number of nodes
        louvainLevel0._phaseI(numRows, isLouvainInit=True) 
        
        for clusterId, cluster in louvainLevel0._clustersLookup.items():
            print()
            self.logger.info("cluserId:{}".format(clusterId))
            self.logger.info(cluster)
            
        expected = {
                    0 : {'clusterId':0, 'numNodes':0, 'weightsInsideCluster':0, 'totalWeight':0},
                    1 : {'clusterId':1, 'numNodes':3, 'weightsInsideCluster':4, 'totalWeight':4},
                    2 : {'clusterId':2, 'numNodes':0, 'weightsInsideCluster':0, 'totalWeight':0}
                    }
            
        self.checkClusters(expected, louvainLevel0._clustersLookup)
            
        self.assertEqual(louvainLevel0._Q, 0.5)
        
        self.logger.info("END\n")
        
    ############################################################
    def testPhaseI(self):
        self.logger.info("BEGIN")
        
        listOfEdges = [(0,1), (1,0),    (0,2), (2,0),    (0,3), (3,0), (1,2), (2,1),
                       (3,4), (4,3) ]
        listOfWeight = [1 for i in listOfEdges]
        louvainLevel0 = Louvain.buildGraph("testPhaseI graph", listOfEdges, listOfWeight)
        louvainLevel0._calculateQ()
        numRows = 5 # number of nodes
        louvainLevel0._phaseI(numRows, isLouvainInit=True)
        
        for clusterId, cluster in louvainLevel0._clustersLookup.items():
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
        
        self.checkClusters(expected, louvainLevel0._clustersLookup)

        c3 = louvainLevel0._clustersLookup[3]
        for n in c3._nodeList:
            self.logger.info("clusterId:{}, nodeId:{}".format(c3._clusterId, n._nodeId))
            
        print('')
        c4 = louvainLevel0._clustersLookup[4]
        for n in c4._nodeList:
            self.logger.info("clusterId:{}, nodeId:{}".format(c4._clusterId, n._nodeId))  
            
        expectedNodesInCluster = {
            3: [0,1,2],
            4: [3,4]
            }
        
        for clusterId in expectedNodesInCluster.keys():
            c = louvainLevel0._clustersLookup[clusterId]  
            nodeIds = [n._nodeId for n in c._nodeList]     
            self.assertEqual(sorted(nodeIds), sorted(expectedNodesInCluster[clusterId]))
        
        self.logger.info("END\n")        
        
    ############################################################        
    def initObjectGraph(self, louvain):
        self.logger.info("BEGIN")
                # init node caches
        for nId in louvain._nodeLookup.keys():
            node = louvain._nodeLookup[nId]
            # because we used _addEdge() instead of addEdges()
            # we need to make sure cache is set up
            node._initKiinCache(louvain._nodeLookup)      
            
        # force nodes to calc cached values
        for nodeId in louvain._nodeLookup.keys():
            node = louvain._nodeLookup[nodeId]
            node.getSumAdjWeights()
            node.getSumOfWeightsInsideCluster(nodeId, louvain._nodeLookup)
            
        # force clusters to calc cached values
        for clusterId in louvain._clustersLookup.keys():
            # run lazy eval
            cluster = louvain._clustersLookup[clusterId]
            cluster.getSumOfWeights()
            cluster.getSumOfWeightsInsideCluster(louvain._nodeLookup) 
            
        self.logger.info("END\n")
        
    ############################################################
    def testPhaseIICreateNewEdges(self):
        self.logger.info("BEGIN")
        listOfEdges = [(0,1), (1,0),    (0,2), (2,0),    (0,3), (3,0), (1,2), (2,1),
                       (3,4), (4,3) ]
        # choose weights to force into two cluster
        listOfWeight = [10,   10,       10,    10,       2,      2,    10,    10 ]
        listOfWeight += [10,  10]
        
        l0 = Louvain.buildGraph("l0", listOfEdges, listOfWeight) 
        
        #  
        # check graph is configured as expected
        #
        self.logger.info("l0 before phase I:\n{}".format(l0)) 
        expectedL0BeforePhaseI = {
            0 : {'clusterId':0 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':22},
            1 : {'clusterId':1 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':20},
            2 : {'clusterId':2 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':20},
            3 : {'clusterId':3 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':12},
            4 : {'clusterId':4 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':10}
            }
        self.checkClusters(expectedL0BeforePhaseI, l0._clustersLookup)        
        
        # check Nodes:
        expectedL0NodesBeforePhaseI = {
            0 : { 'nicd' : {1: {1}, 2: {2}, 3: {3}}, 'wicd': {1: 10, 2: 10, 3: 2} },
            1 : { 'nicd' : {0: {0}, 2: {2}},         'wicd': {0: 10, 2: 10}       },
            2 : { 'nicd' : {0: {0}, 1: {1}},         'wicd': {0: 10, 1: 10}       },
            3 : { 'nicd' : {0: {0}, 4: {4}},         'wicd': {0: 2, 4: 10}        },
            4 : { 'nicd' : {3: {3}},                 'wicd': {3: 10}              }
            }

        for cluster in l0._clustersLookup.values():
            for node in cluster._nodeList:
                self.logger.info("nodeId:{}\n\t _nodesInClusterDict:{}\n\t_weightsInClusterDict:{}"\
                                 .format(node._nodeId, node._nodesInClusterDict, node._weightsInClusterDict))
                expected = expectedL0NodesBeforePhaseI[node._nodeId]
                ret = {'nicd': node._nodesInClusterDict, 'wicd':node._weightsInClusterDict}
                
                self.assertEqual(ret, expected, "clusterId:{} nodeId:{}"\
                                 .format(cluster._clusterId, node._nodeId))
        
        l0._phaseI(numRows=5, isLouvainInit=True)
        self.logger.info("l0 after phase I:\n{}".format(l0)) 
        
        
        AEDWIP        
#         # makes sure is configured as expected
#         expectedL0 = {
#             0 : {'clusterId':0 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':3},
#             1 : {'clusterId':1 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':2},
#             2 : {'clusterId':2 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':2},
#             3 : {'clusterId':3 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':2},
#             4 : {'clusterId':4 ,'numNodes':1 ,'weightsInsideCluster':0 ,'totalWeight':1}
#             }
#         self.checkClusters(expectedL0, l0._clustersLookup)
#         
#         # at this ponit each node is in a separate cluster
#         # create two cluster 
#         c0NodeList = [l0._nodeLookup[0], l0._nodeLookup[1], l0._nodeLookup[2]]
#         c0 = Cluster(0, c0NodeList)  
#         c1NodeList = [l0._nodeLookup[3], l0._nodeLookup[4]]        
#         c1 = Cluster(1, c1NodeList)   
        
        #
        # create a new louvain level that looks like phase I -> phase II -> phase I
        #
        l1 = Louvain("l1", [c0, c1]) 
        self.initObjectGraph(l1)
#         # init node in graph
#         for node in  l1._nodeLookup.values():
#             node._initKiinCache(graphNodesLookup=l1._nodeLookup)
            
        self.logger.info("l1:\n{}".format(l1))
        
        #
        # test edge creation
        #
        l2  = Louvain.buildLouvain('l2', l1)
        l2._phaseIICreateNewEdges()
        
        for nodeId, node in l2._nodeLookup.items():
            self.logger.info("node")
  
        
        self.logger.info("END\n")
                
    ############################################################
    def testPhaseII(self):
        self.logger.info("BEGIN")
        
        listOfEdges = [(0,1), (1,0), (0,2), (2,0), (1,2), (2,1),
                        (0,3), (0,6), (0,7),
                        (3,4), (4,3), (3,5), (5,3),
                        (3,0),
                        (6,7),(7,6), (6,8), (8,6), (6,9), (9,6), (8,9), (9,8), (9,7), (7,9), 
                        (6,0), (7,0)
                       ]
        listOfWeight = [1 for i in listOfEdges]
        louvainLevel0 = Louvain.buildGraph("testPhaseII graph level0", listOfEdges, listOfWeight)
        louvainLevel0._calculateQ()
        numRows = 10 # the number of nodes
        louvainLevel0._phaseI(numRows, isLouvainInit=True)    
        
        expectedAfterPhaseL0_I = {
            0:{'custerId': 0,  'numNodes':0 , 'weightsInsideCluster': 0, 'totalWeight': 0},
            1:{'cluserId': 1,  'numNodes':0 , 'weightsInsideCluster': 0, 'totalWeight': 0},
            2:{'cluserId': 2,  'numNodes':0 , 'weightsInsideCluster': 0, 'totalWeight': 0},
            3:{'cluserId': 3,  'numNodes':0 , 'weightsInsideCluster': 0, 'totalWeight': 0},
            4:{'cluserId': 4,  'numNodes':0 , 'weightsInsideCluster': 0, 'totalWeight': 0},
            5:{'cluserId': 5,  'numNodes':6 , 'weightsInsideCluster':14, 'totalWeight':14},
            6:{'cluserId': 6,  'numNodes':0 , 'weightsInsideCluster': 0, 'totalWeight': 0},
            7:{'cluserId': 7,  'numNodes':0 , 'weightsInsideCluster': 0, 'totalWeight': 0},
            8:{'cluserId': 8,  'numNodes':0 , 'weightsInsideCluster': 0, 'totalWeight': 0},
            9:{'cluserId': 9,  'numNodes':4 , 'weightsInsideCluster':10, 'totalWeight': 12}
            }
        
        self.logger.info("TODO: empty clusters should be pruned, notice")
        
        self.checkClusters(expectedAfterPhaseL0_I, louvainLevel0._clustersLookup)
        
        for cId in [5, 9]:
            nodeIdList = sorted([n._nodeId for n in louvainLevel0._clustersLookup[cId]._nodeList])
            self.logger.info("clusterId:{} nodeList[{}]".format(cId, nodeIdList))
        
        # check phase II
        louvainLevel1 = Louvain.buildLouvain("testPhaseII graph level1", louvainLevel0)
        louvainLevel1._phaseII()
        
        print('')
        self.logger.info("************ check L1 phase II")
        for clusterId, cluster in louvainLevel1._clustersLookup.items():
            self.logger.info("clusterId:{} cluster:{}".format(clusterId,cluster))
            
        expectedAfterPhaseL1_II = {
            5 : {'cluster':5, 'numNodes':1, 'weightsInsideCluster':0, 'totalWeight':2},
            9 : {'cluster':9, 'numNodes':1, 'weightsInsideCluster':0, 'totalWeight':2},
            }
        self.checkClusters(expectedAfterPhaseL1_II, louvainLevel1._clustersLookup)
        
        # check the node caches are set up correctl
        for clusterId, cluster in louvainLevel1._clustersLookup.items():
            for node in cluster._nodeList:
                self.logger.info("clusterId:{} nodeId:{} _weightsInClusterDict:{}"\
                          .format(clusterId, node._nodeId, node._weightsInClusterDict))
                self.logger.info("clusterId:{} nodeId:{} _nodesInClusterDict:{}"\
                          .format(clusterId, node._nodeId, node._nodesInClusterDict))                
        
        # test Louvain algo would run Phase I on louvainLevel2
        # we have to calculate Q before phaseI
        louvainLevel1._calculateQ()
        numRows = 10 # number of nodes        
        louvainLevel1._phaseI(numRows)
        
        print('')
        self.logger.info("************ check L1 after phase I")
        for clusterId, cluster in louvainLevel1._clustersLookup.items():
            self.logger.info("clusterId:{} cluster:{}".format(clusterId,cluster))        
        
        self.logger.info("END\n")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()