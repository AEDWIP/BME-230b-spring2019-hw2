'''
Created on May 16, 2019

@author: andrewdavidson
'''
import logging

############################################################
class Edge(object):
    '''
    TODO
    '''
    logger = logging.getLogger(__name__)
        
    ############################################################        
    def __init__(self, weight, srcId, targetId):
        self._srcId = srcId
        self._targetId = targetId
        self._weight = weight
        
        
    ############################################################                
    def __repr__(self):
        return "srcId:{} targetId:{} weight:{}".format(self._srcId, self._targetId, self._weight)
