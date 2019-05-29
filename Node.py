'''
Created on May 16, 2019

@author: andrewdavidson
'''
import logging

############################################################
class Node(object):
    '''
    public functions:
        __init__(self, clusterId, nodeId)
        __repr__(self)     
        __hash__(self)
        __eq__(self, other)    
        getSumAdjWeights(self)    
        getSumOfWeightsInsideCluster(self, clusterId, graphNodesLookup)
        moveToCluster(self, toClusterId, graphNodesLookup)        
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
        self._adjacentEdgeWeights = 0 #None #0 # edges in and out of cluster
        
        # dynamic programming/caching speed up
        # we are going to need kiin for just about every cluster
        # the value only change when a node is moved to another cluster
        # this does not change Big O, how ever in practice should dramatically improve
        # performance
        self._weightsInClusterDict = {}  # key = clusterId, value = sum of weights in 
        
        # keep a list of all the nodes we are connected to in a given cluster
        # we use this in phase I so we do not test moves that can not improve Q
        self._nodesInClusterDict = {} # key = clusterId, value is a set of nodeIds
        
    ############################################################                
    def __repr__(self):
        str1 = "clusterId:{} nodeId:{} numEdges:{} adjEdgeWeights:{}".format(self._clusterId, 
                                                           self._nodeId, len(self._edgesDict.keys()),
                                                           self._adjacentEdgeWeights)
        str2 = "\n_weightsInClusterDict{}".format(self._weightsInClusterDict)
        str3 = "\n_edgesDict{}".format(self._edgesDict)
        
        return str1 + str2 + str3
    
    ############################################################                    
    def __hash__(self):
        '''
        enable Node objects to be used in sets and dictionary
        '''
        # we assume nodeId is unquie, use this cluster to catch bugs
        return hash((self._clusterId, self._nodeId))

    ############################################################                    
    def __eq__(self, other):
        '''
        enable Node objects to be used in sets and dictionary
        '''        
        if not isinstance(other, type(self)): return NotImplemented
        # out design expects nodeId's to be unique, use cluster to catch bugs
        # checking nodes first is faster
        return self._nodeId == other._nodeId and self._clusterId == other._clusterId
    
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
            self._adjacentEdgeWeights += edge._weight
        else:
            eMsg = "clusterId:{} nodeId:{} edge.targetId:{} was already added".format(self._clusterId, self._nodeId, edge._targetId)
            self.logger.error("ValueError:{}".format(eMsg))
            raise ValueError(eMsg)
        
    ############################################################
    def _initKiinCache(self, graphNodesLookup):
        '''
        enable testing
        '''
        self.logger.debug("BEGIN")
        for key, e in self._edgesDict.items():
            try :
                targetNode = graphNodesLookup[e._targetId]
            except Exception as exc:
                self.logger.warn(exc)
                self.logger.warn("graphNodesLookup[e._targetId] failed: nodeId:{} key:{} edge:{}".format(self._nodeId, key, e))
                continue
            
            targetClusterId = targetNode._clusterId
            self.logger.debug("nodeId:{} e:{} targetNodeClusterId:{}"\
                             .format(self._nodeId, e, targetClusterId))
            if  targetClusterId in self._weightsInClusterDict:
                self._weightsInClusterDict[targetClusterId] += e._weight
                self._nodesInClusterDict[targetClusterId].add(e._targetId)
            else :
                self._weightsInClusterDict[targetClusterId] = e._weight
                self._nodesInClusterDict[targetClusterId] = set()
                self._nodesInClusterDict[targetClusterId].add(e._targetId)

            self.logger.debug("_weightsInClusterDict:{}".format(self._weightsInClusterDict[targetClusterId]))
            self.logger.debug("_nodesInClusterDict:{}".format(self._nodesInClusterDict[targetClusterId]))
                
        self.logger.debug("END\n")
                
        
    ############################################################
    def addEdges(self, listOFEdges,graphNodesLookup):
        '''
        DEPRECATED !!!! you need to construct all the nodes and add all 
        the edges for the entire graph before _initKiinCache
        can raise ValueError
        '''
        for e in listOFEdges:
            self._addEdge(e)

        self._initKiinCache(graphNodesLookup)
          
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
        if self._adjacentEdgeWeights == None:
            self.logger.warning("AEDWIP debug _adjacentEdgeWeight is none!")
            
        return self._adjacentEdgeWeights 
   
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
            
        self.logger.debug("ret:{} clusterId:{}".format(ret, clusterId))
        return ret       
          
    ############################################################
    def _updateAdjNodeMoved(self, adjNodeId, fromClusterId, toClusterId):
        '''
        adjust sufficient stats
        '''
        self.logger.debug("BEGIN")
        
        weight = self.getWeightForEdge(adjNodeId)
        
        # remove the adjNode from our sufficient stats 
        # the adjNode is no longer in the fromCluster
        self._weightsInClusterDict[fromClusterId] -= weight 
        if self._weightsInClusterDict[fromClusterId] == 0 :
            # we are no longer connected to a node in this cluster
            # a tidy obj graph should help track down bug
            check = self._weightsInClusterDict.pop(fromClusterId, None)
            if check == None:
                self.warning("nodeId:{} fromClusterId:{} was missing from _weightsInClusterDict:{}"\
                             .format(self._nodeId, fromClusterId, self._weightsInClusterDict))
                
            check = self._nodesInClusterDict.pop(fromClusterId, None)
            if check == None:
                self.warning("nodeId:{} fromClusterId:{} was missing from _nodesInClusterDict:{}"\
                             .format(self._nodeId, fromClusterId, self._nodesInClusterDict))            
        else:
            self._nodesInClusterDict[fromClusterId].remove(adjNodeId) 
            
        # add the node back to our sufficient stats 
        # it is now in the toCluster
        if not toClusterId in self._weightsInClusterDict:
            self._weightsInClusterDict[toClusterId] = 0
            self._nodesInClusterDict[toClusterId] = set()
            
        self._weightsInClusterDict[toClusterId] += weight
        self._nodesInClusterDict[toClusterId].add(adjNodeId)
        
        self.logger.debug("END\n")
          
    ############################################################
    def moveToCluster(self, toClusterId, graphNodesLookup):
        self.logger.debug("BEGIN nodeId:{} fromClusterId:{} toClusterId:{}"\
                         .format(self._nodeId, self._clusterId, toClusterId))   
        fromClusterId = self._clusterId
        for edge in self._edgesDict.values():
            adjNodeId = edge._targetId
            adjNode = graphNodesLookup[adjNodeId]
            adjNode._updateAdjNodeMoved(self._nodeId, fromClusterId, toClusterId)
            
            
        # update our selvs
        # the nodes we are connected to did not change
        self._clusterId = toClusterId
        
        
        self.logger.debug("END\n")  
        
