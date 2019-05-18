'''
Created on May 16, 2019

@author: andrewdavidson
'''
import logging

############################################################
class Node(object):
    '''
    classdocs
    '''

    logger = logging.getLogger(__name__)

    ############################################################
    def __init__(self, clusterId, nodeId):
        '''
        Constructor
        '''
        
        self._clusterId = clusterId
        self._nodeId = nodeId
        self._edges = {} # key is the edge's target id, value is the edge obj
        self._adjcentEdgeWeights = 0 # edges in and out of cluster
        self._weightsInCluster = None # edges only between nodes in cluster

        
    ############################################################                
    def __repr__(self):
        return "clusterId:{} nodeId:{} numEdges:{} adjEdgeWeights:{}".format(self._clusterId, 
                                                           self._nodeId, len(self._edges.keys()),
                                                           self._adjcentEdgeWeights)
    
    ############################################################
    def addEdge(self, edge):
        '''
        can raise ValueError
        '''
        
        if not edge._targetId in self._edges:
            self._edges[edge._targetId] = edge
            self._adjcentEdgeWeights += edge._weight
        else:
            eMsg = "clusterId:{} nodeId:{} edge.targetId:{} was already added".format(self._clusterId, self._nodeId, edge._targetId)
            self.logger.error("ValueError:{}".format(eMsg))
            raise ValueError(eMsg)
          
    ############################################################
    def _getEdges(self):
        return [v for k,v in self._edges.items()]
        
    ############################################################
    def getSumAdjWeights(self):
        '''
        This is the 'Ki' term in Louvain paper 
        "Fast unfolding of communities in large networks"
        '''
        return self._adjcentEdgeWeights

    ############################################################
    def getSumOfWeightsInsideCluster(self, graphNodesLookup):
        '''
        This is the 'Kin' term in Louvain paper 
        "Fast unfolding of communities in large networks"
        
        arguments:
            graphNodesLookup:
                a dictionary or set of nodes in graph. keys should be 
                node ids. 
        '''
        if not self._weightsInCluster:
            self._weightsInCluster  = 0
            for e in self._edges:
                targetId = e._targetId
                if targetId in graphNodesLookup:
                    n = graphNodesLookup[e._targetId]
                    if self._clusterId == n._clusterId:
                        self._weightsInCluster += e._weight
                else:
                    eMsg = "targetId:{} not found in graphNodesLookup".format(targetId)
                    self.logger.error(eMsg)
                    raise ValueError(eMsg)
            
        return self._weightsInCluster        
    
    ############################################################
    def getM(self):   
        '''
        The nodes contribution to m in 
        "Fast unfolding of communities in large networks"
        
        returns 1/2 * self.getSumAdjWeights()
        '''
        return 0.5 * self.getSumAdjWeights()
         
   
    ############################################################
    def getWeightForEdge(self, edgeTargetId):
        '''
        This is the 'Aij' term in Louvain paper 
        "Fast unfolding of communities in large networks"
        
        can raise ValueError
        '''
        ret = 0
        if edgeTargetId in self._edges:
            ret = self._edges[edgeTargetId]._weight
        else:
            eMsg = "clusterId:{} nodeId:{} edge.targetId:{} missing".format(self._clusterId, self._nodeId, edgeTargetId)
            self.logger.error("ValueError:{}".format(eMsg))
            raise ValueError(eMsg)        
        
        return ret
    

        