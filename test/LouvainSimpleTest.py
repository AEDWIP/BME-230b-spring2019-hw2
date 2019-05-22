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
from scanpy.tools._louvain import louvain


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

        # create second cluster graph
        n3 = Node(clusterId="c1", nodeId=3)
        e6 = Edge(weight=1.0, srcId=3, targetId=4)
        n3._addEdge(e6)

        n4 = Node(clusterId="c1", nodeId=4)
        e6 = Edge(weight=1.0, srcId=4, targetId=3)
        n4._addEdge(e6)
        
        cluster1 = Cluster(clusterId="c1", nodeList=[n3, n4])
        self.assertEqual(1, cluster1._getM())
        clusters = [cluster0, cluster1]
        
        # you can not move a node to a cluster if the node is not
        # connected to something in the cluster
        # there would not gain in Q
        # create and edge between a node in c0 and c2

        ea = Edge(weight=1.0, srcId=n0._nodeId, targetId=n3._nodeId)
        eb = Edge(weight=1.0, srcId=n3._nodeId, targetId=n0._nodeId)
        n0._addEdge(ea)
        n3._addEdge(eb)        

        level0 = Louvain("simple", [cluster0, cluster1] )  
              
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
        level0 = Louvain("testNode", [cluster0, cluster1])
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
        level0 = Louvain("testQ", clusterList)
        self.assertEqual(5, level0._getM())
        
        self.logger.info("level0._Q:{}".format(level0._Q))
        self.assertAlmostEqual(level0._Q, 0.44)        
        self.logger.info("END\n")
        
    ############################################################
    def testChangeInQSlow(self):
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
        
        louvain1 = Louvain("changeInQ1", [cluster1, cluster2])
        
        # calculate modularity of original graph
        self.logger.info("louvain1._Q:{}".format(louvain1._Q))
        self.assertEqual(louvain1._Q, 0.5599999999999999)
        
        # move node 2 from cluster 2 to cluster 1
        n2._clusterId = "c1"
        cluster1 = Cluster(clusterId="1", nodeList=[n0, n1, n2, n3])
        cluster2 = Cluster(clusterId="2", nodeList=[n4, n5])
        
        # calculate modularity
        louvain2 = Louvain("changeInQ2", [cluster1, cluster2])
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
            self.assertEqual(c._totalWeight, expectedData["sigmaTotal"], eMsg)
            self.assertEqual(c._weightsInsideCluster, expectedData["sigmaIn"], eMsg)
        self.logger.info("END\n")    
    

    ############################################################    
    def checkKiinStats(self, msg, clusters, expectedKiinData):
        self.logger.info("BEGIN")  
        self.logger.info(msg) 
        for c in clusters:
            cid = c._clusterId
            expectedData = expectedKiinData[cid]
            for n in c._nodeList:
                nid = n._nodeId
                key = "k{}in".format(nid)
                eMsg = "clusterId:{} nodeId:{} key:{}".format(cid, nid, key)
                self.assertEqual(n._weightsInClusterDict[cid], expectedData[key], eMsg)
        
        self.logger.info("END\n")  

    ############################################################    
    def checkNodeStats(self, msg, clusters, expectedNodeData):
        self.logger.info("BEGIN") 
        self.logger.info(msg)   
        for c in clusters:
            cid = c._clusterId
            for n in c._nodeList:
                nid = n._nodeId
                eMsg = "clusterId:{} nodeId:{}".format(cid, nid)
                expectedData = expectedNodeData[nid]
                self.assertEqual(len(n._edgesDict.keys()), expectedData["numEdges"], eMsg)
                self.assertEqual(n._clusterId, expectedData["clusterId"], eMsg)
                self.assertEqual(n._adjcentEdgeWeights, expectedData["adjEdgeWeights"], eMsg)

        self.logger.info("END\n")    
    
    ############################################################        
    def testMove(self):
        self.logger.info("BEGIN")    
               
        graph = self.createSimpleGraph()
        louvainLevel = graph[0]
        clusters = graph[1]
        nodes = graph[2]
        graphNodesLookup = graph[4]
        
#         # you can not move a node to a cluster if the node is not
#         # connected to something in the cluster
#         # there would not gain in Q
#         # create and edge between a node in c0 and c2
#         na = nodes[0]
#         nb = nodes[3]
#         self.assertNotEqual(na._clusterId, nb._clusterId)
#         ea = Edge(weight=1.0, srcId=na._nodeId, targetId=nb._nodeId)
#         eb = Edge(weight=1.0, srcId=nb._nodeId, targetId=na._nodeId)
#         na._addEdge(ea)
#         nb._addEdge(eb)
        
        # make sure the graph is set up the way we expect
        self.logger.info("********* make sure graph is configured correctly")
        for n in nodes:
            print()
            self.logger.info("nodeId:{} clusterId:{}".format(n._nodeId, n._clusterId))
            for eid, edge in n._edgesDict.items():
                self.logger.info(edge)
        
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
        beforeExpectedClusterData = {
            "c0": {"numNodes":3, 'sigmaIn':6, 'sigmaTotal':7},
            "c1": {"numNodes":2, 'sigmaIn':2, 'sigmaTotal':3}
            }   
        self.checkClusterStats("********** before move", clusters, beforeExpectedClusterData)

        beforeExpectedKiinData = {
            "c0": {"k0in":2, "k1in":2, "k2in":2, "k3in":1, "k4in":0},
            "c1": {"k0in":1, "k1in":0, "k2in":0, "k3in":1, "k4in":1}            
            }
        self.checkKiinStats("********** before move", clusters, beforeExpectedKiinData)
        
        beforeMoveQ = louvainLevel._Q
        self.logger.info("beforeMoveQ:{}".format(beforeMoveQ))
           
        # test move
        c0 = clusters[0]
        c1 = clusters[1]
        n0 = nodes[0]
        c0.moveNode(c1, n0, graphNodesLookup)
        
        # check nodes
        self.logger.info("")
        afterExpectedNodeData = {
                    0:{'clusterId':'c1', 'numEdges':3, 'adjEdgeWeights':3.0},
                    1:{'clusterId':'c0', 'numEdges':2, 'adjEdgeWeights':2.0},
                    2:{'clusterId':'c0', 'numEdges':2, 'adjEdgeWeights':2.0},
                    3:{'clusterId':'c1', 'numEdges':2, 'adjEdgeWeights':2.0},
                    4:{'clusterId':'c1', 'numEdges':1, 'adjEdgeWeights':1.0}}
        self.checkNodeStats("********** after move", clusters, afterExpectedNodeData)    
      
#         for n in nodes:
#             self.logger.info("after node:{} \n_weightsInClusterDict:\n{}\n".format(n, n._weightsInClusterDict))
        
        # check kiin
        afterExpectedKiinData = {
            "c0": {"k0in":2, "k1in":1, "k2in":1, "k3in":0, "k4in":0},
            "c1": {"k0in":2, "k1in":1, "k2in":1, "k3in":2, "k4in":1}            
            }
        self.checkKiinStats("********** after move", clusters, afterExpectedKiinData)    
        
        # check cluster     
        for c in clusters:
            self.logger.info(c) 
              
        afterExpectedClusterData = {
            "c0": {"numNodes":2, 'sigmaIn':2, 'sigmaTotal':4},
            "c1": {"numNodes":3, 'sigmaIn':4, 'sigmaTotal':6}
            }  
        self.checkClusterStats("********** after move", clusters, afterExpectedClusterData)
 
        # what is Q after move?
        louvainLevel._calculateQ()
        afterMoveQ = louvainLevel._Q
        changeInQ = afterMoveQ - beforeMoveQ
        self.logger.info("changeInQ:{} afterMoveQ:{} before:{}".format(changeInQ, afterMoveQ, beforeMoveQ))
        expectedBeforeMoveQ = 0.44
        expecedAfterMoveQ = 0.36
        self.assertAlmostEqual(afterMoveQ, expecedAfterMoveQ)
        self.assertAlmostEqual(beforeMoveQ, expectedBeforeMoveQ)

        
        # expectedChangeInQ assumes our implementation of _calculateQ() and moveNode()
        # are correct. I did not verify this value by hand
        expectedChangeInQ = -0.08  #-0.28125  
        self.assertAlmostEqual(changeInQ, expectedChangeInQ)
        
        self.logger.warn("TODO: AEDWIP: try and calculate what the change is by fast formula")
        
        self.logger.info("END\n")            


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
