'''
Created on May 23, 2019

@author: andrewdavidson
'''
from Cluster import Cluster
from Edge import Edge
import logging
from louvain import Louvain
from Node import Node
import numpy as np
import pandas as pd
from setupLogging import setupLogging
import unittest
import scanpy.api as sc
from euclid_knn import KnnG
from timeit import default_timer as timer
from datetime import timedelta

class LouvainScanpyIntegrationTest(unittest.TestCase):

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
    def testGetClusterAssigments(self):
        self.logger.info("BEGIN")
        
        # build lovain level tree
        listOfEdges = [(0,1), (1,0), (0,2), (2,0), (1,2), (2,1),
                        (0,3), (0,6), (0,7),
                        (3,4), (4,3), (3,5), (5,3),
                        (3,0),
                        (6,7),(7,6), (6,8), (8,6), (6,9), (9,6), (8,9), (9,8), (9,7), (7,9), 
                        (6,0), (7,0)
                       ]
        listOfWeight = [1 for i in listOfEdges]
        
        # build out init graph and calculate Q
        louvainLevel0 = Louvain.buildGraph("level0", listOfEdges, listOfWeight)
        louvainLevel0._calculateQ() 
               
        # run phase I: find best cluster assignments
        louvainLevel0._phaseI(isLouvainInit=True) 
        
        ll_0_ClusterAssigments = louvainLevel0.getClusterAssigments()
        self.logger.info("level0 cluster assignments:\n{}".format(ll_0_ClusterAssigments))  
                
        # create next level and run phaseII
        # phase II consolidates clusters found in previous level
        louvainLevel1 = Louvain.buildLouvain("level1", louvainLevel0)
        louvainLevel1._phaseII()
        louvainLevel1._calculateQ()
        
        # lets assume this is the top level.
        louvainLevel1._phaseI()     
 
        ll_1_ClusterAssigments = louvainLevel1.getClusterAssigments()
        self.logger.info("level1 cluster assignments:\n{}".format(ll_1_ClusterAssigments))  

        self.logger.info("**************** check for side effects output should be same as above")
        ll_0_ClusterAssigments = louvainLevel0.getClusterAssigments()
        self.assertEqual(ll_0_ClusterAssigments, {5: [0, 1, 2, 3, 4, 5], 9: [9, 6, 7, 8]})

        ll_1_ClusterAssigments = louvainLevel1.getClusterAssigments()
        self.assertEqual(ll_1_ClusterAssigments, {9: [9, 6, 7, 8, 0, 1, 2, 3, 4, 5]})

        self.logger.info("END\n")
           
    ############################################################
    def testRun(self):
        self.logger.info("BEGIN")
        
        # build lovain level tree
        listOfEdges = [(0,1), (1,0), (0,2), (2,0), (1,2), (2,1),
                        (0,3), (0,6), (0,7),
                        (3,4), (4,3), (3,5), (5,3),
                        (3,0),
                        (6,7),(7,6), (6,8), (8,6), (6,9), (9,6), (8,9), (9,8), (9,7), (7,9), 
                        (6,0), (7,0)
                       ]
        listOfWeight = [1 for i in listOfEdges]
        root = Louvain.run(listOfEdges, listOfWeight)
        rootClusterAssigments = root.getClusterAssigments()
        self.logger.info("louvainId: {} root cluster assigments:\n{}".format(root._louvainId, rootClusterAssigments))
        self.assertEqual(rootClusterAssigments, {9: [9, 6, 7, 8, 0, 1, 2, 3, 4, 5]})

    ############################################################
    def testTransformClusterAssigments(self):
        '''
        convert cluster assignments to format scanpy expects
        '''
        self.logger.info("BEGIN")
        rootAssigment = {5: [0, 1, 2, 3, 4, 5], 9: [9, 6, 7, 8]}
        flatTuples = []
        for clusterId, nodeList in rootAssigment.items():
            for nodeId in nodeList:
                flatTuples.append( (nodeId, clusterId) )
        self.logger.info("flatTuples:{}".format(flatTuples))
                
        # sort by node id
        sortedTuples = sorted(flatTuples,  key=lambda x: x[0])
        self.logger.info("sortedTuples:{}".format(sortedTuples))
        
        clusterAssigments = [t[1] for t in sortedTuples]
        self.logger.info("clusterAssigments:{}".format(clusterAssigments))
        self.assertEqual(clusterAssigments, [5, 5, 5, 5, 5, 5, 9, 9, 9, 9])
        
        idx = [chr(65 + i) for i in range(10)]
        assignmentPS = pd.Series(data=clusterAssigments, index=idx)
        self.logger.info("assignmentPS:\n{}".format(assignmentPS))
        
        self.logger.info("END\n")
        
    ############################################################
    @unittest.skip("skip this test it takes over 3 mins")
    def testAdata(self):
        '''
        The real deal
        '''
        self.logger.info("BEGIN") 
        anndata = sc.read("../PBMC.merged.h5ad")
               
        xxx = anndata.uns['neighbors']['connectivities']
        
        
# anndata.uns['neighbors']['connectivities'] csr_matrix
# shape <class 'tuple'>: (15476, 15476)
# 
# adata.obs['louvain']
# Series: index
# data_3p-AAACCTGAGCATCATC-0     9
# data_3p-AAACCTGAGCTAGTGG-0     5
# <class 'tuple'>: (15476,)

        #knn takes about 3 or 4 min
        # run our implementation of nearest neighboors and update anndata
# todo try running with out knng maybe adata has values already save time
        KnnG(anndata, n_neighbors=12, runPCA=True, nPC=50)
        
        self.logger.info("begin Louvain.runWithAdata")
        start = timer()
        root = Louvain.runWithAdata(anndata)
        end = timer()
        self.logger.info("Louvain.runWithAdata execution time:{}"\
                         .format(timedelta(seconds=end-start)))        
        
        
        self.logger.info("END\n")        
        
    ############################################################
    def testGetAdataInput(self):
        '''
        TODO
        '''
        self.logger.info("BEGIN")     
        anndata = sc.read("../PBMC.merged.h5ad")
        
        # run our implementation of nearest neighboors and update anndata
        KnnG(anndata, n_neighbors=12, runPCA=True, nPC=50)        
               
        adjacency = anndata.uns['neighbors']['connectivities']
        sources, targets = adjacency.nonzero()
        listOfWeights = adjacency[sources, targets]
        if isinstance(listOfWeights, np.matrix):
            # Return listOfWeights as a flattened np.array
            # https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.matrix.A1.html
            listOfWeights = listOfWeights.A1
            
        listOfEdges = list(zip(sources, targets))
  
        self.logger.info("END\n")  
               
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()