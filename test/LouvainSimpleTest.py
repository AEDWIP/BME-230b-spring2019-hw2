'''
Created on May 17, 2019

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
class LouvainSimpleTest(unittest.TestCase):
    setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)

    ############################################################
    def setUp(self):
        pass

    ############################################################
    def tearDown(self):
        # make sure all the logs are flushed
        # if assert we will not see partial test log entries
        logging.shutdown()
        
    ############################################################
    def createSimpleGraph(self):
        '''
        creates two disjoint graphs, one is a triangle, the other is a pair of nodes
        connected by a single edge
        
        creates two cluster. one for each of the disjoint graphs
        
        all edge weights are 1
        
        returns (level0Louvain, 
               [cluster0, cluster1], 
               [n0, n1, n2, n3, n4], 
               [e0, e1, e2, e3, e4, e5, e6])
        
        '''
        self.logger.info("BEGIN")
        n0 = Node(clusterId="c0", nodeId=0)
        n1 = Node(clusterId="c0", nodeId=1)
        n2 = Node(clusterId="c0", nodeId=2)

        # undirected  triangle graph
        e0 = Edge(weight=1.0, srcId=0, targetId=1)
        n0._addEdge(e0)

        e1 = Edge(weight=1.0, srcId=1, targetId=0)
        n1._addEdge(e1)
                
        e2 = Edge(weight=1.0, srcId=0, targetId=2)
        n0._addEdge(e2)   
        
        e3 = Edge(weight=1.0, srcId=2, targetId=0) 
        n2._addEdge(e3)    
        
        e4 = Edge(weight=1.0, srcId=1, targetId=2)
        n1._addEdge(e4)
        
        e5 = Edge(weight=1.0, srcId=2, targetId=1) 
        n2._addEdge(e5) 
        
        # create  cluster0
        cluster0 = Cluster(clusterId="c0", nodeList=[n0, n1, n2])

        # create disjoint graph
        n3 = Node(clusterId="c1", nodeId=3)
        e6 = Edge(weight=1.0, srcId=3, targetId=4)
        n3._addEdge(e6)

        n4 = Node(clusterId="c1", nodeId=4)
        e6 = Edge(weight=1.0, srcId=4, targetId=3)
        n4._addEdge(e6)
        
        cluster1 = Cluster(clusterId="c1", nodeList=[n3, n4])
        self.assertEqual(1, cluster1._getM())
        clusters = [cluster0, cluster1]

        level0 = Louvain([cluster0, cluster1] )  
              
        self.logger.info("END\n")
        
        nodeList = [n0, n1, n2, n3, n4]
        graphNodesLookup = { n._nodeId:n for n in nodeList}
        
#         for n in nodeList:
#             # because we used _addEdge() instead of addEdges()
#             # we need to make sure cache is set up
#             n._initKiinCache(graphNodesLookup)
            
        ret = (level0, 
               clusters, 
               nodeList, 
               [e0, e1, e2, e3, e4, e5, e6], 
               graphNodesLookup)
        
        return (ret)

    ############################################################
    def testNode(self):
        self.logger.info("BEGIN")
        
        n0 = Node(clusterId="c0", nodeId=0)
        n1 = Node(clusterId="c0", nodeId=1)
        n2 = Node(clusterId="c0", nodeId=2)

        # undirected  triangle graph
        e0 = Edge(weight=1.0, srcId=0, targetId=1)
        n0._addEdge(e0)

        e1 = Edge(weight=1.0, srcId=1, targetId=0)
        n1._addEdge(e1)
        
        self.assertEqual(1, n0.getSumAdjWeights())
        self.assertEqual(1, n1.getSumAdjWeights())
        
        e2 = Edge(weight=1.0, srcId=0, targetId=2)
        n0._addEdge(e2)   
        
        e3 = Edge(weight=1.0, srcId=2, targetId=0) 
        n2._addEdge(e3)    
        
        # test print functions
        self.logger.info("e3:{}".format(e3))        
        
        self.assertEqual(2, n0.getSumAdjWeights())
        self.assertEqual(1, n2.getSumAdjWeights())
        
        # test print functions
        self.logger.info("n2:{}".format(n2))

        e4 = Edge(weight=1.0, srcId=1, targetId=2)
        n1._addEdge(e4)
        
        e5 = Edge(weight=1.0, srcId=2, targetId=1) 
        n2._addEdge(e5) 
        
        self.assertEqual(2, n1.getSumAdjWeights())
        self.assertEqual(2, n2.getSumAdjWeights())        

        # create  cluster0
        cluster0 = Cluster(clusterId="c0", nodeList=[n0, n1, n2])
        self.assertEqual(3, cluster0._getM())

        # test print functions
        self.logger.info("cluster0:{}".format(cluster0))        

        # create disjoint graph
        n3 = Node(clusterId="c1", nodeId=3)
        e6 = Edge(weight=1.0, srcId=3, targetId=4)
        n3._addEdge(e6)

        n4 = Node(clusterId="c1", nodeId=4)
        e6 = Edge(weight=1.0, srcId=4, targetId=3)
        n4._addEdge(e6)
        
        cluster1 = Cluster(clusterId="c1", nodeList=[n3, n4])
        self.assertEqual(1, cluster1._getM())

        # test modularity calculation
        level0 = Louvain([cluster0, cluster1])
        self.assertEqual(4, level0._getM())
        
        self.logger.info("level0._Q:{}".format(level0._Q))
        self.assertEqual(level0._Q, 0.59375)
        
        # test 
        
        self.logger.info("END\n")

    ############################################################
    def testQ(self):
        self.logger.info("BEGIN")
        
        # test modularity calculation
        ret = self.createSimpleGraph()
        clusterList = ret[1]
        self.logger.info("clusterList:\n{}".format(clusterList))
        level0 = Louvain(clusterList)
        self.assertEqual(4, level0._getM())
        
        self.logger.info("level0._Q:{}".format(level0._Q))
        self.assertEqual(level0._Q, 0.59375)        
        self.logger.info("END\n")
        
    ############################################################
    def testChangeInQ(self):
        self.logger.info("BEGIN")
        
        n0 = Node(clusterId="c1", nodeId=0)
        n1 = Node(clusterId="c1", nodeId=1)
        n3 = Node(clusterId="c1", nodeId=3)
        
        e1 = Edge(weight=1.0, srcId=0, targetId=1) 
        n0._addEdge(e1)
        e2 = Edge(weight=1.0, srcId=0, targetId=2) 
        n0._addEdge(e2)
        e3 = Edge(weight=1.0, srcId=0, targetId=3) 
        n0._addEdge(e3)

        e4 = Edge(weight=1.0, srcId=1, targetId=0) 
        n1._addEdge(e4)
        

        e5 = Edge(weight=1.0, srcId=3, targetId=0) 
        n3._addEdge(e5)
        
        cluster1 = Cluster(clusterId="1", nodeList=[n0, n1, n3])

        n2 = Node(clusterId="c2", nodeId=2)
        e6 = Edge(weight=1.0, srcId=2, targetId=0) 
        n2._addEdge(e6)

        
        n4 = Node(clusterId="c2", nodeId=4)
        n5 = Node(clusterId="c2", nodeId=5)
        e7 = Edge(weight=1.0, srcId=4, targetId=5) 
        n4._addEdge(e7)
        e8 = Edge(weight=1.0, srcId=5, targetId=4) 
        n5._addEdge(e8)
        
        e9 = Edge(weight=1.0, srcId=4, targetId=2) 
        n4._addEdge(e9)
        
        e10 = Edge(weight=1.0, srcId=2, targetId=4) 
        n2._addEdge(e10)
        
        cluster2 = Cluster(clusterId="2", nodeList=[n2, n4, n5])
        
        louvain1 = Louvain([cluster1, cluster2])
        
        # calculate modularity of original graph
        self.logger.info("louvain1._Q:{}".format(louvain1._Q))
        self.assertEqual(louvain1._Q, 0.5599999999999999)
        
        # move node 2 from cluster 2 to cluster 1
        n2._clusterId = "c1"
        cluster1 = Cluster(clusterId="1", nodeList=[n0, n1, n2, n3])
        cluster2 = Cluster(clusterId="2", nodeList=[n4, n5])
        
        # calculate modularity
        louvain2 = Louvain([cluster1, cluster2])
        self.logger.info("louvain2._Q:{}".format(louvain2._Q))
        self.assertEqual(louvain2._Q, 0.5199999999999999)

        self.logger.info("change in modularity:{}".format(louvain1._Q - louvain2._Q))
    
    ############################################################
    def checkClusterStats(self, msg, clusters, expectedClusterData):
        self.logger.info("BEGIN")  
        self.logger.info(msg)  
        for c in clusters:
            eMsg = "clusterId:{}".format(c._clusterId)
            expectedData = expectedClusterData[c._clusterId]
            self.assertEqual(len(c._nodeList), expectedData["numNodes"], eMsg)
            self.assertEqual(c._weightsInsideCluster, expectedData["sigmaIn"], eMsg)
            self.assertEqual(c._totalWeight, expectedData["sigmaTotal"], eMsg)

        self.logger.info("END\n")    
    
    ############################################################
    def testMove(self):
        self.logger.info("BEGIN")    
               
        graph = self.createSimpleGraph()
        clusters = graph[1]
        nodes = graph[2]
        graphNodesLookup = graph[4]
        
        # you can not move a node to a cluster if the node is not
        # connected to something in the cluster
        # there would not gain in Q
        # create and edge between a node in c0 and c2
        na = nodes[0]
        nb = nodes[3]
        self.assertNotEqual(na._clusterId, nb._clusterId)
        ea = Edge(weight=1.0, srcId=na._nodeId, targetId=nb._nodeId)
        eb = Edge(weight=1.0, srcId=nb._nodeId, targetId=na._nodeId)
        na._addEdge(ea)
        nb._addEdge(eb)
        
        for n in nodes:
            # because we used _addEdge() instead of addEdges()
            # we need to make sure cache is set up
            n._initKiinCache(graphNodesLookup)        
            
        self.logger.info("")                    
        for n in nodes:
            # for lazy evaluation to run
            n.getSumAdjWeights()
            n.getSumOfWeightsInsideCluster(n._clusterId, graphNodesLookup)
            self.logger.info("node:{}".format(n))
            
        self.logger.info("")            
        for c in clusters:
            # run lazy eval
            c.getSumOfWeights()
            c.getSumOfWeightsInsideCluster(graphNodesLookup)
            self.logger.info("cluster:{}".format(c))
               
        # test if cluster is init correctly
        beforeExpectedClusterData = np.array([
                                            [],
                                            [],
                                             ])
        beforeExpectedClusterData = {
            "c0": {"numNodes":3, 'sigmaIn':6, 'sigmaTotal':7},
            "c1": {"numNodes":2, 'sigmaIn':2, 'sigmaTotal':3}
            }   
        
        self.checkClusterStats("before move", clusters, beforeExpectedClusterData)
        
# #         self.assertEqual(c0._weightsInsideCluster, 6.0)
# #         self.assertEqual(c0._totalWeight, 7.0)
# #         
# #         c1 = clusters[1]
# #         self.assertEqual(c1._weightsInsideCluster, 2.0)
# #         self.assertEqual(c1._totalWeight, 3.0)        
#         
# #         n3 = nodes[3]
# #         self.logger.info("before move 3 is buggy node:{} \n_weightsInClusterDict:\n{}\n".format(n3, n3._weightsInClusterDict))
#         
#         # is data buggy?
#         for n in nodes:
#             self.logger.info("before move buggy? node:{} \n_weightsInClusterDict:\n{}\n".format(n, n._weightsInClusterDict))
#             

        # test move
        c0 = clusters[0]
        c1 = clusters[1]
        c0.moveNode(c1, na, graphNodesLookup)
        
        afterExpectedClusterData = {
            "c0": {"numNodes":2, 'sigmaIn':2, 'sigmaTotal':4},
            "c1": {"numNodes":3, 'sigmaIn':4, 'sigmaTotal':6}
            }  
        
#         self.logger.info("after move:c0:{}".format(c0))
#         self.logger.info("after move:c1:{}".format(c1))
#         
        # check cluster
        self.checkClusterStats("after move", clusters, afterExpectedClusterData)
        
#         self.assertEqual(c0._totalWeight, 4.0)        
#         self.assertEqual(c1._totalWeight, 6.0)  
                
        # TODO 
#         self.assertEqual(c0._weightsInsideCluster, 2.0)
#         self.assertEqual(c1._weightsInsideCluster, 4.0)
        
        # check nodes
        for n in nodes:
            self.logger.info("node:{} \n_weightsInClusterDict:\n{}\n".format(n, n._weightsInClusterDict))
            
#         self.assertEqual(nodes[0]._weightsInClusterDict, {'c0': 2.0, 'c1': 1.0})
#         self.assertEqual(nodes[1]._weightsInClusterDict, {'c0': 1.0})
#         self.assertEqual(nodes[2]._weightsInClusterDict, {'c0': 1.0})
#         self.assertEqual(nodes[3]._weightsInClusterDict, {'c1': 2.0, 'c0': 0.0})
# 
#         self.assertEqual(nodes[4]._weightsInClusterDict, {'c1': 1.0})     
        
        self.logger.info("END\n")            


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()