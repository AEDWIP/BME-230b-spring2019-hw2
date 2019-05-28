'''
Created on May 28, 2019

@author: andrewdavidson
'''
# from Cluster import Cluster
# from Edge import Edge
import logging
from louvain import Louvain
# from Node import Node
# import numpy as np
from setupLogging import setupLogging
import unittest


class TestAEDWIP(unittest.TestCase):

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
    def testAEDWIP(self):
        self.logger.info("BEGIN")
        
        # full connected graph
        # with three super nodes
        # should converge to 3 clusters
        listOfEdges   = [ (0,i) for i in range(1,5)]
        listOfWeights = [1 for i in range(1,5)]
        
        listOfEdges   += [(5, i) for i in range(6,10)]
        listOfWeights += [1 for i in range (6,10)]
        
        listOfEdges   += [(10, i) for i in range(11,15)]
        listOfWeights += [1 for i in range (11,15)]        
        numRows = 15
        
        # add edges connecting the super nodes
        listOfEdges +=  [(0,5), (5,0),  (5,10), (10, 5)]
        listOfWeights += [1,     1,     1,       1]
            
        self.logger.info(listOfEdges)
        self.logger.info(listOfWeights)
        
        self.logger.info("len(listOfEdges):{} len(listOfWeights):{}"\
                         .format(len(listOfEdges), len(listOfWeights)))
        
        louvainLevel0 = Louvain.buildGraph("level 0", listOfEdges, listOfWeights)
        self.logger.info("louvainLevel0 after buildGraph():\n{}".format(louvainLevel0))
         
        self.logger.info("*********** did we insert a blank line" )        
        self.logger.info('')    
        self.logger.info("*********** did we insert a blank line" )        

        # run phase I
        louvainLevel0._phaseI(numRows, isLouvainInit=True) # TODO: can probably get rid of isLouvainInit
        self.logger.info("after phase I() louvainLevel0:\n{}".format(louvainLevel0))
        l0Assignments = louvainLevel0.getClusterAssigments()
        self.logger.info("l0Assigments cluster assignments:\n{}".format(l0Assignments))        
          
#         # check Q
#         louvainLevel0._calculateQ()
#         self.logger.info("after phase I   louvainLevel0._Q:{}".format(louvainLevel0._Q))
#            
#         # build next level
#         louvainLevel1 = Louvain.buildLouvain("level 1 ", louvainLevel0)
#         self.logger.info("after buildLouvain louvainLevel1\n{}".format(louvainLevel1))
#          
#         # phase II
#         louvainLevel1._phaseII(isLouvainInit=False) # TODO: can probably get rid of isLouvainInit)
#         self.logger.info("after phaseII() louvainLevel1  this log line looks funnny:\n{}".format(louvainLevel1)) 
#         l1Assignments = louvainLevel1.getClusterAssigments()
#         self.logger.info("louvainLevel1 cluster assignments:\n{}".format(l1Assignments))          

        self.logger.info("END\n")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()