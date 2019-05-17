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
    def __init__(self, targetId, weight):
        self._targetId = targetId
        self._weight = weight
        
