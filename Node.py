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
        self._edgesDict = {} # key is the edge's target id, value is the edge obj
        self._adjcentEdgeWeights = 0 # edges in and out of cluster
        
        # dynamic programming/caching speed up
        # we are going to need kiin for just about every cluster
        # the value only change when a node is moved to another cluster
        # this does not change Big O, how ever in practice should dramatically improve
        # performance
        self._weightsInClusterDict = {} # key = clusterId, value = sum of weights in 

        
    ############################################################                
    def __repr__(self):
        return "clusterId:{} nodeId:{} numEdges:{} adjEdgeWeights:{}".format(self._clusterId, 
                                                           self._nodeId, len(self._edgesDict.keys()),
                                                           self._adjcentEdgeWeights)
    
    ############################################################
    def _addEdge(self, edge):
        '''
        use addEdges(listofEdges) in production
        
        -addEdge can only be used for simple unit tests. 
        you need to make sure you call _initKiinCache after the last 
        - addEdge() else move() will not work correct
        
        can raise ValueError
        '''
        
        if not edge._targetId in self._edgesDict:
            self._edgesDict[edge._targetId] = edge
            self._adjcentEdgeWeights += edge._weight
        else:
            eMsg = "clusterId:{} nodeId:{} edge.targetId:{} was already added".format(self._clusterId, self._nodeId, edge._targetId)
            self.logger.error("ValueError:{}".format(eMsg))
            raise ValueError(eMsg)
        
    ############################################################
    def _initKiinCache(self, graphNodesLookup):
        '''
        enable testing
        '''
        self.logger.info("BEGIN")
        for key, e in self._edgesDict.items():
            targetNode = graphNodesLookup[e._targetId]
            targetClusterId = targetNode._clusterId
            self.logger.info("nodeId:{} e:{} targetNodeClusterId:{}"\
                             .format(self._nodeId, e, targetClusterId))
            if  targetClusterId in self._weightsInClusterDict:
                self._weightsInClusterDict[targetClusterId] += e._weight
            else :
                self._weightsInClusterDict[targetClusterId] = e._weight
                
            self.logger.info("_weightsInClusterDict[{}]:{}".format(targetClusterId, self._weightsInClusterDict[targetClusterId]))
                
        self.logger.info("END\n")
                
        
    ############################################################
    def addEdges(self, listOFEdges):
        '''
        can raise ValueError
        '''
        for e in listOFEdges:
            self._addEdge(e)
#             
#             # initialize kiin cache
#             if e._targetId in self._weightsInClusterDict:
#                 self._weightsInClusterDict[e._targetId] += e._weight
#             else :
#                 self._weightsInClusterDict[e._targetId] = e._weight

        self._initKiinCache()
          
    ############################################################
    def _getEdges(self):
        return [v for k,v in self._edgesDict.items()]
        
    ############################################################
    def getSumAdjWeights(self):
        '''
        This is the 'Ki' term in Louvain paper 
        "Fast unfolding of communities in large networks"
        '''
        # this only changes if an edge is added
        return self._adjcentEdgeWeights 
    
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
        if edgeTargetId in self._edgesDict:
            ret = self._edgesDict[edgeTargetId]._weight
        else:
            eMsg = "clusterId:{} nodeId:{} edge.targetId:{} missing".format(self._clusterId, self._nodeId, edgeTargetId)
            self.logger.error("ValueError:{}".format(eMsg))
            raise ValueError(eMsg)        
        
        return ret
 
    ############################################################
    def getSumOfWeightsInsideCluster(self, clusterId, graphNodesLookup):
        '''
        This is the 'Kin' term in Louvain paper 
        "Fast unfolding of communities in large networks"
        
        arguments:
            graphNodesLookup:
                a dictionary or set of nodes in graph. keys should be 
                node ids. 
        '''
        
        ret = 0
        if clusterId in self._weightsInClusterDict:
            ret = self._weightsInClusterDict[clusterId]
        else :
            eMsg = "nodeId:{} clusterId:{} not found in _weightsInClusterDict:\n{}"\
                .format(self._nodeId, clusterId, self._weightsInClusterDict)
            self.logger.error(eMsg)
            raise ValueError(eMsg)
        
#         else :
#             self._weightsInClusterDict  = 0
#             ret = 0
#             for e in self._edgesDict:
#                 targetId = e._targetId
#                 if targetId in graphNodesLookup:
#                     n = graphNodesLookup[e._targetId]
#                     if self._clusterId == n._clusterId:
#                         ret += e._weight
#                 else:
#                     eMsg = "targetId:{} not found in graphNodesLookup".format(targetId)
#                     self.logger.error(eMsg)
#                     raise ValueError(eMsg)
#                 
#             self._weightsInClusterDict[clusterId] = ret
            
        self.logger.info("ret:{} clusterId:{}".format(ret, clusterId))
        return ret       
    
#      ############################################################
#     def _addedToCluster(self, clusterId, graphNodesLookup):
#         '''
#         use move()!!!!!!
#         add() and remove() must be called in a particular order        
#         '''
#         
#         
#     ############################################################
#     def _removedFromCluster(self, graphNodesLookup):
#         '''
#         use move()!!!
#         add() and remove() must be called in a particular order
#         '''

        
    ############################################################
    def moveToCluster(self, clusterId, graphNodesLookup):
        '''
        TODO:
        '''
        # move does not effect self kiin for fromCluster or toCluster        
        for e in self._edgesDict:
            n = graphNodesLookup[e._targetId]
            n._adjustKiin(nodeId=self._nodeId, 
                         fromCluster=self._clusterId,
                         toCluster=clusterId,
                         weight=e._weight)
            
        self._clusterId = clusterId
        
    ############################################################
    def _adjustKiin(self, nodeId, fromClusterId, toClusterId, weight):
        '''
        a node we are connected to was moved into a new cluster.
        '''
        needToDecriment = fromClusterId in self._weightsInClusterDict
        needToIncrement = toClusterId in self._weightsInClusterDict
        
        if not (needToDecriment or needToIncrement):
            eMsg = "self.nodeId{} in not connected to nodeId:{} fromClusterId:{} toClusterId:{}"\
                .format(self._nodeId, nodeId, fromClusterId, toClusterId)
            self.logger.error(eMsg)
            raise ValueError(eMsg)
        
        if needToDecriment:
            self._weightsInClusterDict[fromClusterId] -= weight
            
        if needToIncrement:
            self._weightsInClusterDict[fromClusterId] += weight
            
        