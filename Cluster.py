'''
Created on May 16, 2019

@author: andrewdavidson
'''

import logging

############################################################
class Cluster(object):
    '''
    public functions:
        __init__(self, clusterId, nodeList)
        __repr__(self)
        getSumOfWeightsInsideCluster(self, graphNodesLookup)
        getSumOfWeights(self)
        moveNode(self, targetCluster, node, graphNodesLookup, isLouvainInit=False)    
    '''

    logger = logging.getLogger(__name__)
    
    ############################################################      
    def __init__(self, clusterId, nodeList):
        self._clusterId = clusterId
        self._nodeList = nodeList
        self._weightsInsideCluster = None # 0 # 'sigma in'
        self._totalWeight = None # 0 # 'sigma tot'
        self.getSumOfWeights() 
        
    ############################################################                
    def __repr__(self):
        return "clusterId:{} numNodes:{} :weightsInsideCluster:{} totalWeight:{}".\
            format(self._clusterId, len(self._nodeList),self._weightsInsideCluster, self._totalWeight)
    
    ############################################################
    def _getEdges(self):
        ret = []
        for n in self._nodeList:
            ret += n._getEdges()
            
        self.logger.debug("c:{} ret:\n{}".format(self._clusterId, ret))
        return ret
        
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
                
        TODO: AEDWIP: FIXME: if all the nodes wind up in a single cluster the value should be 1/2
        we do not use SigmaIn so deal with this as time permits
        '''
        
        if not self._weightsInsideCluster:
            self._weightsInsideCluster  = 0
            for n in self._nodeList:
                self.logger.debug("clusterId:{} nodeId:{}".format(self._clusterId, n._nodeId))
                kiin = n.getSumOfWeightsInsideCluster(self._clusterId, graphNodesLookup)
                self._weightsInsideCluster += kiin
            
        return self._weightsInsideCluster                   

    ############################################################
    def getSumOfWeights(self):
        '''
        This is the 'Sigma total' term in Louvain paper 
        "Fast unfolding of communities in large networks"
        
        TODO: AEDWIP: FIXME: this value is not correct. If all nodes
        get move into a single cluster the value should be 1/2
        we do not use sigmaTotal os not an issue fix as time permits
        '''
        if not self._totalWeight:
            self._totalWeight = 0
            for n in self._nodeList:
                self._totalWeight += n.getSumAdjWeights()
            
        return self._totalWeight
        
    ############################################################
    def _removeNode(self, node, graphNodesLookup, isLouvainInit=False):
        '''
        in production use move()
        
        arguments:
            isLouvainInit: 
                the initialization step puts each node in a separate cluster. by definition
                node._weightsInClusterDict[self._clusterId] is un defined.
                if not the initialization step, this should be an error 
                
                default value is False
        '''
        
        self._totalWeight -= node.getSumAdjWeights()   
        kiin = node.getSumOfWeightsInsideCluster(self._clusterId, graphNodesLookup)
        self.logger.debug("clusterId:{} nodeId:{}  _weightsInsideCluster:{} kiin:{}"\
                .format(self._clusterId, node._nodeId, self._weightsInsideCluster, kiin))  
                      

        if not isLouvainInit:
            # account for edges that get transformed from inside our cluster to
            # between our cluster. remember we model links between nodes a pair
            # of directed edges
            
            # there is not requirement that we have links to other nodes
            # in our cluster
            # TODO: get rid of isLouvainInit it is not a special case
            if self._clusterId in node._weightsInClusterDict:
                w = 2 * node._weightsInClusterDict[self._clusterId]
                self._weightsInsideCluster -= w
        
        self._nodeList.remove(node)  
          
    ############################################################
    def _addNode(self, node, targetClusterId, graphNodesLookup):
        '''
        TODO
        
        in production use move()

        '''
        self.logger.debug("node._adjacentEdgeWeights:{}".format(node._adjacentEdgeWeights))
        if not node:
            self.loggger.warning("AEDWIP DEBUG node is none!!")
            
        if self._totalWeight == None:
            self.loggger.warning("AEDWIP DEBUG _totalWeight is none!!")
            
        self._totalWeight += node.getSumAdjWeights()
        kiin = node.getSumOfWeightsInsideCluster(targetClusterId, graphNodesLookup)
        self.logger.debug("clusterId:{} nodeId:{} targetClusterId:{} _weightsInsideCluster:{} kiin:{}"\
                .format(self._clusterId, node._nodeId, targetClusterId, self._weightsInsideCluster, kiin))
        
        # we gain twice. there is an edge between the node being moved
        # and a node in the target cluster. There is also a node from the cluster with and 
        # edge back to the node being moved
        self._weightsInsideCluster += 2 * kiin
        
        self._nodeList.append(node)
        
    ############################################################
    def moveNode(self, targetCluster, node, graphNodesLookup, isLouvainInit=False):
        '''
        TODO
        
        arguments:
            isLouvainInit: 
                the initialization step puts each node in a separate cluster. by definition
                node._weightsInClusterDict[self._clusterId] is undefined.
                if not the initialization step, this should be an error 
                
                default value is False        
        '''
        self.logger.debug("clusterId:{} nodeId:{} targetClusterId:{}"\
                         .format(self._clusterId, node._nodeId, targetCluster._clusterId))
        
        targetClusterId = targetCluster._clusterId        
        targetCluster._addNode(node, targetClusterId, graphNodesLookup)
        self._removeNode(node, graphNodesLookup, isLouvainInit)    
        node.moveToCluster(targetCluster._clusterId, graphNodesLookup )
            
