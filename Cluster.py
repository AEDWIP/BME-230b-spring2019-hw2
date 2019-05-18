'''
Created on May 16, 2019

@author: andrewdavidson
'''

import logging

############################################################
class Cluster(object):
    '''
    TODO:
    '''

    logger = logging.getLogger(__name__)
    
    ############################################################      
    def __init__(self, clusterId, nodeList):
        self._clusterId = clusterId
        self._nodeList = nodeList
        self._weightsInsideCluster = None
        self._totalWeight = None
        
    ############################################################                
    def __repr__(self):
        return "clusterId:{} numNodes:{} :weightsInsideCluster:{}".format(self.clusterId, len(self._nodeList),
                                                 self._weightsInsideCluster )
    
    ############################################################
    def _getEdges(self):
        ret = []
        for n in self._nodeList:
            ret += n._getEdges()
            
        self.logger.info("c:{} ret:\n{}".format(self._clusterId, ret))
        return ret
    
    ############################################################
    def _getM(self):
        '''
        the partial m term in the Louvain paper
        "Fast unfolding of communities in large networks"
        
        returns 1/2 the sum of all edges in the the cluster
        '''
        m = 0
        for node in self._nodeList:
            m += node.getM() 
            
        return m       
    
    ############################################################
    def _getNodes(self):
        '''
        returns a list of nodes
        '''
        return self._nodeList

    ############################################################
    def getSumOfWeightsInsideCluster(self, graphNodesLookup):
        '''
        This is the 'Sigma in term' in Louvain paper 
        "Fast unfolding of communities in large networks"
        
        arguments:
            graphNodesLookup:
                a dictionary or set of nodes in graph. keys should be 
                node ids. 
        '''
        if not self._weightsInCluster:
            self._weightsInCluster  = 0
            for n in self._nodeList:
                self._weightsInsideCluster += n.getSumOfWeightsInsideCluster(graphNodesLookup)
            
        return self._weightsInCluster                    

    ############################################################
    def getSumOfWeights(self):
        '''
        This is the 'Sigma total' term in Louvain paper 
        "Fast unfolding of communities in large networks"
        '''
        if not self._totalWeight:
            self._totalWeight = 0
            for n in self._nodeList:
                self._totalWeight += n.getSumAdjWeights()
            
        return self.self._totalWeight
        