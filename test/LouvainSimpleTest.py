'''
Created on May 17, 2019

@author: andrewdavidson
'''
import logging
from setupLogging import setupLogging
import unittest
from Edge import Edge
from Node import Node
from Cluster import Cluster
from louvain import Louvain

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
        n0.addEdge(e0)

        e1 = Edge(weight=1.0, srcId=1, targetId=0)
        n1.addEdge(e1)
                
        e2 = Edge(weight=1.0, srcId=0, targetId=2)
        n0.addEdge(e2)   
        
        e3 = Edge(weight=1.0, srcId=2, targetId=0) 
        n2.addEdge(e3)    
        
        e4 = Edge(weight=1.0, srcId=1, targetId=2)
        n1.addEdge(e4)
        
        e5 = Edge(weight=1.0, srcId=2, targetId=1) 
        n2.addEdge(e5) 
        
        # create  cluster0
        cluster0 = Cluster(clusterId="c0", nodeList=[n0, n1, n2])

        # create disjoint graph
        n3 = Node(clusterId="c1", nodeId=3)
        e6 = Edge(weight=1.0, srcId=3, targetId=4)
        n3.addEdge(e6)

        n4 = Node(clusterId="c1", nodeId=4)
        e6 = Edge(weight=1.0, srcId=4, targetId=3)
        n4.addEdge(e6)
        
        cluster1 = Cluster(clusterId="c1", nodeList=[n3, n4])
        self.assertEqual(1, cluster1._getM())

        level0 = Louvain([cluster0, cluster1] )  
              
        self.logger.info("END\n")
        
        ret = (level0, 
               [cluster0, cluster1], 
               [n0, n1, n2, n3, n4], 
               [e0, e1, e2, e3, e4, e5, e6])
        
        return (ret)


    ############################################################
    def testNode(self):
        self.logger.info("BEGIN")
        
        n0 = Node(clusterId="c0", nodeId=0)
        n1 = Node(clusterId="c0", nodeId=1)
        n2 = Node(clusterId="c0", nodeId=2)

        # undirected  triangle graph
        e0 = Edge(weight=1.0, srcId=0, targetId=1)
        n0.addEdge(e0)

        e1 = Edge(weight=1.0, srcId=1, targetId=0)
        n1.addEdge(e1)
        
        self.assertEqual(1, n0.getSumAdjWeights())
        self.assertEqual(1, n1.getSumAdjWeights())
        
        e2 = Edge(weight=1.0, srcId=0, targetId=2)
        n0.addEdge(e2)   
        
        e3 = Edge(weight=1.0, srcId=2, targetId=0) 
        n2.addEdge(e3)    
        
        # test print functions
        self.logger.info("e3:{}".format(e3))        
        
        self.assertEqual(2, n0.getSumAdjWeights())
        self.assertEqual(1, n2.getSumAdjWeights())
        
        # test print functions
        self.logger.info("n2:{}".format(n2))

        e4 = Edge(weight=1.0, srcId=1, targetId=2)
        n1.addEdge(e4)
        
        e5 = Edge(weight=1.0, srcId=2, targetId=1) 
        n2.addEdge(e5) 
        
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
        n3.addEdge(e6)

        n4 = Node(clusterId="c1", nodeId=4)
        e6 = Edge(weight=1.0, srcId=4, targetId=3)
        n4.addEdge(e6)
        
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
        n0.addEdge(e1)
        e2 = Edge(weight=1.0, srcId=0, targetId=2) 
        n0.addEdge(e2)
        e3 = Edge(weight=1.0, srcId=0, targetId=3) 
        n0.addEdge(e3)

        e4 = Edge(weight=1.0, srcId=1, targetId=0) 
        n1.addEdge(e4)
        

        e5 = Edge(weight=1.0, srcId=3, targetId=0) 
        n3.addEdge(e5)
        
        cluster1 = Cluster(clusterId="1", nodeList=[n0, n1, n3])

        n2 = Node(clusterId="c2", nodeId=2)
        e6 = Edge(weight=1.0, srcId=2, targetId=0) 
        n2.addEdge(e6)

        
        n4 = Node(clusterId="c2", nodeId=4)
        n5 = Node(clusterId="c2", nodeId=5)
        e7 = Edge(weight=1.0, srcId=4, targetId=5) 
        n4.addEdge(e7)
        e8 = Edge(weight=1.0, srcId=5, targetId=4) 
        n5.addEdge(e8)
        
        e9 = Edge(weight=1.0, srcId=4, targetId=2) 
        n4.addEdge(e9)
        
        e10 = Edge(weight=1.0, srcId=2, targetId=4) 
        n2.addEdge(e10)
        
        cluster2 = Cluster(clusterId="2", nodeList=[n2, n4, n5])
        
        louvain = Louvain([cluster1, cluster2])
        
        # Q:0.47000000000000003
        self.logger.info("Q:{}".format(louvain._Q))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()