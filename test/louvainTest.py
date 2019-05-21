'''
Created on May 16, 2019

@author: andrewdavidson
'''
import igraph as ig
import logging
from louvain import Louvain
import numpy as np
# import scanpy.api as sc
# from scipy.spatial.distance import pdist
from setupLogging import setupLogging
import unittest

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
        louvainLevel0 = Louvain.buildGraph(listOfEdges, listOfWeight)
        self.logger.info("Q:{}".format(louvainLevel0._Q))
        self.logger.info("END\n")
        
    ############################################################
    def testChangeInModulartiy(self):
        self.logger.info("BEGIN")
        
        listOfEdges = [(0,1), (0,2), (0,3), (2,4), (4,5) ]
        listOfWeight = [1 for i in listOfEdges]
        
        self.logger.info("END\n")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()