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
        self._adjcentEdgeWeight = 0
        
    ############################################################                
    def __repr__(self):
        return "clusterId:{} nodeId:{} numEdges:{} adjEdgeWeight:{}".format(self._clusterId, 
                                                           self._nodeId, len(self._edges.keys()),
                                                           self._adjcentEdgeWeight)
    
    ############################################################
    def addEdge(self, edge):
        '''
        can raise ValueError
        '''
        
        if not edge._targetId in self._edges:
            self._edges[edge._targetId] = edge
            self._adjcentEdgeWeight += edge._weight
        else:
            eMsg = "clusterId:{} nodeId:{} edge.targetId:{} was already added".format(self._clusterId, self._nodeId, edge._targetId)
            self.logger.error("ValueError:{}".format(eMsg))
            raise ValueError(eMsg)
          
    ############################################################
    def _getEdges(self):
        return [v for k,v in self._edges.items()]
        
    ############################################################
    def getSumAdjWeight(self):
        '''
        This is the Ki term in Louvain paper 
        "Fast unfolding of communities in large networks"
        '''
        return self._adjcentEdgeWeight

    
    ############################################################
    def getM(self):   
        '''
        The nodes contribution to m in 
        "Fast unfolding of communities in large networks"
        
        returns 1/2 * self.getSumAdjWeight()
        '''
        return 0.5 * self.getSumAdjWeight()
         
   
    ############################################################
    def getWeightForEdge(self, edgeTargetId):
        '''
        This is the Aij term in Louvain paper 
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
    
   
        