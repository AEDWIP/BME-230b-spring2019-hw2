'''
Created on May 16, 2019

@author: andrewdavidson
'''
import igraph as ig
import logging
import numpy as np
import scanpy.api as sc
from scipy.spatial.distance import pdist
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
    def createSimpleGraph(self):
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
    def testCreateSimpleGraph(self):
        g = self.createSimpleGraph()
        self.logger.info("g:\n:{}".format(g))
        self.logger.info("g.es['weight']:\n:{}".format(g.es['weight']))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()