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

############################################################
class Test(unittest.TestCase):
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
    def testNode(self):
        self.logger.info("BEGIN")
        
        n0 = Node(clusterId="c0", nodeId=0)
        n1 = Node(clusterId="c0", nodeId=1)
        n2 = Node(clusterId="c0", nodeId=2)

        # undirected  triangle graph
        e0 = Edge(weight=1, srcId=0, targetId=1)
        n0.addEdge(e0)

        e1 = Edge(weight=1, srcId=1, targetId=0)
        n1.addEdge(e1)
        
        self.assertEqual(1, n0.getSumAdjWeight())
        self.assertEqual(1, n1.getSumAdjWeight())
        
        e2 = Edge(weight=1, srcId=0, targetId=2)
        n0.addEdge(e2)   
        
        e3 = Edge(weight=1, srcId=2, targetId=0) 
        n2.addEdge(e3)    
        
        self.assertEqual(2, n0.getSumAdjWeight())
        self.assertEqual(1, n2.getSumAdjWeight())

        e4 = Edge(weight=1, srcId=2, targetId=1)
        n1.addEdge(e4)
        
        e5 = Edge(weight=1, srcId=1, targetId=2) 
        n2.addEdge(e5) 
        
        self.assertEqual(2, n1.getSumAdjWeight())
        self.assertEqual(2, n2.getSumAdjWeight())        

        #
        # create  cluster
        #
        cluster0 = Cluster(clusterId="c0", nodeList=[n0, n1, n2])
        self.assertEqual(3, cluster0._getM())
#         


        
        self.logger.info("END\n")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()