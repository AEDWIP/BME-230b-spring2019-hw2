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
        self._edges = []
        self._adjEdgeWeight = None
        
    ############################################################
    def addEdge(self, edge):
        self._edges.append(edge)
          
    ############################################################
    def sumAdjWeight(self):
        ret = self._adjEdgeWeight
        if not self._adjEdgeWeight:
            w = 0
            for edge in self._edges :
                w += edge._weight
                
            self._adjEdgeWeight = w
                
            
        
        