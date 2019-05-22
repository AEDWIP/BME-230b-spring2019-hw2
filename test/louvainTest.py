'''
Created on May 16, 2019

@author: andrewdavidson
'''
import igraph as ig
import logging
from louvain import Louvain
import LouvainSimpleTest
import numpy as np
# import scanpy.api as sc
# from scipy.spatial.distance import pdist
from setupLogging import setupLogging
import unittest
from Edge import Edge
from Node import Node
from Cluster import Cluster

############################################################
class LouvainTest(unittest.TestCase):
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
    def createSimpleIGraph(self):
        '''
        example from iGraph tutorial 
        https://igraph.org/python/doc/tutorial/tutorial.html
        '''
        g = ig.Graph([(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)])
        # add attributes to the vertices
        g.vs["name"] = \
           ["Alice", "Bob", "Claire", "Dennis", "Esther", "Frank", "George"]
        g.vs["age"] = [25, 31, 18, 47, 22, 23, 50]
        g.vs["gender"] = ["f", "m", "f", "m", "f", "m", "m"]
        
        # add attributes to the edges
        g.es["is_formal"] = [False, False, True, True, True, False, True, False, False]
        
        # add weights
        # to keep things simple assume edges have a weight of one
        # with anndata the weights will be distances between cell values
        g.es['weight'] = weights = [1 for i in range(len(g.vs)) ]
        g.es['weight'] = weights
        
        return g
        
    ############################################################
    def testCreateSimpleIGraph(self):
        self.logger.info("BEGIN")
        g = self.createSimpleIGraph()
        self.logger.info("g:\n:{}".format(g))
        self.logger.info("g.es['weight']:\n:{}".format(g.es['weight']))
        self.logger.info("END\n")

    ############################################################
    def testBootStrapModularity(self):
        self.logger.info("BEGIN")
        #g.community_multilevel()
        
        g = self.createSimpleIGraph()
        
        # modularity membership is a list with length = number of nodes
        # the value in the list corresponds to the cluster the node is
        ml = [i for i in range(g.vcount())]
        self.logger.info("membership:{}".format(ml))
        expectedModularity = g.modularity(ml)
        self.logger.info("iGraph Modularity:{}".format(expectedModularity))
        self.logger.warn("iGraph Modularity can not be used to test bootstrap")
        self.logger.warn("the cluster only have a single node. no edge is inside how come modularity is not 0")
        
        # test out code
        listOfEdges = [ e.tuple for e in g.es]
        self.logger.info("listOfEdges:\n{}".format(listOfEdges))
        
        listOfWeight = list(g.es['weight']) #[e['weights'] for e in g.es]
        self.logger.info("listOfWeight:\n{}".format(listOfWeight))

#         louvainLevel0 = Louvain(None)
#         louvainLevel0.bootStrapInit(listOfEdges, listOfWeight)
        louvainLevel0 = Louvain.buildGraph("testBootStrapModularity", listOfEdges, listOfWeight)
        self.logger.info("Q:{}".format(louvainLevel0._Q))
        self.logger.info("END\n")
        
    ############################################################    
    def createChangeQGraph(self):
        # create cluster 0
        n0 = Node(clusterId="c0", nodeId=0)
        n1 = Node(clusterId="c0", nodeId=1)
        n2 = Node(clusterId="c0", nodeId=2)        
        n3 = Node(clusterId="c1", nodeId=3)
        n4 = Node(clusterId="c1", nodeId=4)

        # 0 - 1
        e0 = Edge(weight=1.0, srcId=0, targetId=1)
        n0._addEdge(e0)
        e1 = Edge(weight=1.0, srcId=1, targetId=0)
        n1._addEdge(e1)
               
        # 0 - 2 
        e2 = Edge(weight=1.0, srcId=0, targetId=2)
        n0._addEdge(e2)   
        e3 = Edge(weight=1.0, srcId=2, targetId=0) 
        n2._addEdge(e3)    
        
        # 0 - 3
        # edge between clusters
        e4 = Edge(weight=1.0, srcId=0, targetId=3)
        n0._addEdge(e4)
        e5 = Edge(weight=1.0, srcId=3, targetId=0) 
        n3._addEdge(e5) 
        
        cluster0 = Cluster(clusterId="c0", nodeList=[n0, n1, n2])

        # creat cluster 1
        #n5 = Node(clusterId="c1", nodeId=5)

        # 4 - 3
        e6 = Edge(weight=1.0, srcId=4, targetId=3)
        n4._addEdge(e6)
        e7 = Edge(weight=1.0, srcId=3, targetId=4)
        n3._addEdge(e7)
        
        cluster1 = Cluster(clusterId="c1", nodeList=[n3, n4])
        clusters = [cluster0, cluster1]
        
        e8 = Edge(weight=1.0, srcId=2, targetId=1)
        n2._addEdge(e8)
        e9 = Edge(weight=1.0, srcId=1, targetId=2)
        n1._addEdge(e9)
        
        edgeList = [e0, e1, e2, e3, e4, e5, e6, e7, e8, e9]
        i = 1
        for e in edgeList:
            if i % 2 :
                print()
                
            i += 1
            self.logger.info(e)  
        print()      
        
        nodeList = [n0, n1, n2, n3, n4]
        graphNodesLookup = { n._nodeId:n for n in nodeList}
        
        for n in nodeList:
            # because we used _addEdge() instead of addEdges()
            # we need to make sure cache is set up
            n._initKiinCache(graphNodesLookup)        
            
        self.logger.debug("")                    
        for n in nodeList:
            # for lazy evaluation to run
            n.getSumAdjWeights()
            n.getSumOfWeightsInsideCluster(n._clusterId, graphNodesLookup)
            self.logger.debug("node:{}".format(n))
            
        self.logger.debug("")            
        for c in clusters:
            # run lazy eval
            c.getSumOfWeights()
            c.getSumOfWeightsInsideCluster(graphNodesLookup)
            self.logger.info("cluster:{}".format(c))        
            
        level0 = Louvain("changeQGraph", [cluster0, cluster1] )  
        ret = (level0, 
               clusters, 
               nodeList, 
               edgeList, 
               graphNodesLookup)
        
        self.logger.info("END\n")        
        return (ret)
           
    ############################################################
    def testChangeInModulartiy(self):
        self.logger.info("BEGIN")
        
#         this does not build graph the way we need for this test
#         listOfEdges = [(0,1), (0,2), (0,3), (2,4), (4,5) ]
#         listOfWeight = [1 for i in listOfEdges]
#         louvainLevel0 = Louvain.buildGraph(listOfEdges, listOfWeight)

        graph = self.createChangeQGraph()
        louvain = graph[0]
        clusters = graph[1]
        nodesList = graph[2]
        graphNodesLookup = graph[4]
        self.logger.info("louvainLevel0:{}".format(louvain))
        
        # check if nodeSets are correct
        for n in nodesList:
            print()
            self.logger.info("nodeId:{}".format(n._nodeId))
            for clusterId in n._nodesInClusterDict.keys():
                self.logger.info("\t clusterId: {} set:{}".format(clusterId, n._nodesInClusterDict[clusterId]))
        
        # check modularity before move
        self.assertAlmostEqual(louvain._Q, 0.44)
        
        n0 = nodesList[0]
        fromCluster = clusters[0]
        targetCluster = clusters[1]
        
        # test change if node was removed from cluster
        removeChange = louvain.changeInModularityIfNodeRemoved(n0, fromCluster)
        self.logger.info("removeChange:{}".format(removeChange))
        self.assertAlmostEqual(removeChange, 0.16)
        
        # test what change would be if we moved n2 from cluster 0 to cluster 1
        ret =louvain.modularityGainIfMove(fromCluster, targetCluster, n0)
        
        expectedChangeInQ = -0.08
        self.logger.info("modularityGainIfMove:{} expected:{}".format(ret, expectedChangeInQ))
        self.assertEqual(ret, expectedChangeInQ)

        self.logger.info("END\n")
        
    ############################################################
    def testChangeInModulartiyTrivalGraph(self):
        self.logger.info("BEGIN")
        
        # create a trival graph.
        listOfEdges = [(0,1), (1,2)]
        listOfWeight = [1 for i in listOfEdges]
        louvainLevel0 = Louvain.buildGraph("trival graph", listOfEdges, listOfWeight)

        self.logger.info("louvainLevel0:{}".format(louvainLevel0))
        
        # make sure graph is set u as expected
        for nid,node in louvainLevel0._nodeLookup.items():
            print()
            for eid, edge in node._edgesDict.items():
                self.logger.info(edge)
        print()
        
        # check modularity before move
        beforeMoveQ = louvainLevel0._Q
        self.assertEqual(beforeMoveQ, 0.0) 

        n1 = louvainLevel0._nodeLookup[1]
        fromCluster = louvainLevel0._clusters[1]
        targetCluster = louvainLevel0._clusters[2]
         
        predictedChangeInQ =louvainLevel0.modularityGainIfMove(fromCluster, targetCluster, n1)
        self.logger.info("predicted changeInQ:{}".format(predictedChangeInQ))
        
        # move
        fromCluster.moveNode(targetCluster, n1, louvainLevel0._nodeLookup, isLouvainInit=True)
        
        # calculate Q
        louvainLevel0._calculateQ()
        afterMoveQ = louvainLevel0._Q
        expectedChangeInQ = afterMoveQ - beforeMoveQ
        self.logger.info("expectedChangeInQ:{} afterMoveQ:{} before:{}"\
                         .format(expectedChangeInQ, afterMoveQ, beforeMoveQ))
        self.logger.info("predicted changeInQ:{}".format(predictedChangeInQ))

#         
#         expectedChangeInQ = -0.04
#         self.logger.info("modularityGainIfMove:{} expected:{}".format(ret, expectedChangeInQ))
#         self.assertEqual(ret, expectedChangeInQ)

        self.logger.info("END\n")        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()