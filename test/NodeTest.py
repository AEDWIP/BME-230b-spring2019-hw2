'''
Created on May 17, 2019

@author: andrewdavidson
'''
import logging
from setupLogging import setupLogging
import unittest
from Edge import Edge
from Node import Node

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

        # undirected graph
        e0 = Edge(weight=1, srcId=0, targetId=1)
        e1 = Edge(weight=1, srcId=1, targetId=0)
        n0.addEdge(e0)
        n1.addEdge(e1)
        
        self.assertEqual(1, n0.getSumAdjWeight())
        self.assertEqual(1, n1.getSumAdjWeight())

#         
#         e1 = Edge(weight=1, srcId=0, targetId=2)
#         e2 = Edge(weight=1, srcId=0, targetId=1)
#         e1 = Edge(weight=1, srcId=0, targetId=2)
#         
#         # in cartesian coordinates e3 weight should be > 1
#         e2 = Edge(weight=1, srcId=1, targetId=2)
#         


        
        self.logger.info("END\n")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()