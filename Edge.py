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
    
    ############################################################                    
    def __hash__(self):
        '''
        enable Edge objects to be used in sets and dictionary
        '''
#         uniqueId = str(self._weight) + str(self._srcId) + str(self._targetId)
        return hash((self._weight, self._srcId, self._targetId))

    ############################################################                    
    def __eq__(self, other):
        '''
        enable Edge objects to be used in sets and dictionary
        '''        
        if not isinstance(other, type(self)): return NotImplemented
        return self._weight == other._weight and self._targetId == other._targetId and self._srcId == other._srcId    
