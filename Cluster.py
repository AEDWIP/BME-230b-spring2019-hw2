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
        self._weightsInsideCluster = None # 'sigma in'
        self._totalWeight = None # 'sigma tot'
        
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
        
        if not self._weightsInsideCluster:
            self._weightsInsideCluster  = 0
            for n in self._nodeList:
                self.logger.info("clusterId:{} nodeId:{}".format(self._clusterId, n._nodeId))
                kiin = n.getSumOfWeightsInsideCluster(self._clusterId, graphNodesLookup)
                if not kiin: # TODO: AEDWIP
                    self.logger.info("kiin WTF?")
#                 elif not self._weightsInsideCluster: value of zero drops us into this block
#                     self.logger.info("_weightsInsideCluster WTF? if value is zero okay:{}".format(self._weightsInsideCluster))

                self._weightsInsideCluster += kiin
            
        return self._weightsInsideCluster                   

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
            
        return self._totalWeight
        
    ############################################################
    def _removeNode(self, node, graphNodesLookup):
        '''
        in production use move()
        '''

        self._totalWeight -= node.getSumAdjWeights()   
        kiin = node.getSumOfWeightsInsideCluster(self._clusterId, graphNodesLookup)
        self._weightsInsideCluster -= kiin
        
        # TODO maintain an index remove in linear, index is constant time
        self._nodeList.remove(node)    
    ############################################################
    def _addNode(self, node, targetClusterId, graphNodesLookup):
        '''
        TODO
        
        in production use move()

        '''
        if self._clusterId == targetClusterId:
            self.logger.info("do not pass target id it confusing") # this is the case
        else:
            self.logger.info("we need the target id")
            
        self._totalWeight += node.getSumAdjWeights()
        kiin = node.getSumOfWeightsInsideCluster(targetClusterId, graphNodesLookup)
        self._weightsInsideCluster += kiin
        
        self._nodeList.append(node)
        
    ############################################################
    def moveNode(self, targetCluster, node, graphNodesLookup):
        '''
        TODO
        '''
        targetClusterId = targetCluster._clusterId        
        targetCluster._addNode(node, targetClusterId, graphNodesLookup)
        
        self._removeNode(node, graphNodesLookup)

        node.moveToCluster(targetCluster._clusterId, graphNodesLookup)