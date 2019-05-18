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
        
    ############################################################                
    def __repr__(self):
        return "clusterId:{} numNodes:{}".format(self.clusterId, len(self._nodeList))
    
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
            
