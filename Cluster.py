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
        return "clusterId:{} numNodes:{} :weightsInsideCluster:{} totalWeight:{}".\
            format(self._clusterId, len(self._nodeList),self._weightsInsideCluster, self._totalWeight)
    
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
        This is the 'Sigma in' term: in Louvain paper 
        "Fast unfolding of communities in large networks"
        
        arguments:
            graphNodesLookup:
                a dictionary or set of nodes in graph. keys should be 
                node ids. 
        '''
        if not self._weightsInCluster:
            self._weightsInCluster  = 0
            for n in self._nodeList:
                kiin = n.getSumOfWeightsInsideCluster(self._clusterId, graphNodesLookup)
                self._weightsInsideCluster += kiin
            
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
        
    ############################################################
    def _removeNode(self, node):
        '''
        TODO
        '''
        recalculate self getSumOfWeights
        recalculate self _weightsInCluster 

        node.removedFromCluster(self._clusterId)
        
    ############################################################
    def _addNode(self, node):
        '''
        TODO
        '''
        recalculate self getSumOfWeights   
        recalculate self _weightsInCluster 
        node.addedToCluster(self._clusterId) 
        

