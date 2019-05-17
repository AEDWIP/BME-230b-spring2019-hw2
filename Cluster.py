'''
Created on May 16, 2019

@author: andrewdavidson
'''
import igraph as ig
import igraph
from posix import access
import logging

get rid of igraph and just roll your own
############################################################
class Cluster(object):
    '''
    classdocs
    '''

    logger = logging.getLogger(__name__)
    
    ############################################################      
    def __init__(self, nodeList):
        self._nodeList = nodeList
        
    ############################################################
    def _calculateM(self):
        '''
        the partial m term in the Louvain paper
        "Fast unfolding of communities in large networks"
        
        returns 1/2 the sum of all edges in the the cluster
        '''
        m = 0
        for node in self._nodeList:
            m += node._calculateM() 
            
        return m       
        
#        '''
#         arguments
#             graph: 
#                 ref: https://igraph.org/python/doc/igraph.Graph-class.html
#             nodeList: 
#                 a subset of nodes in graph g that are to be assigned to the this cluster
#                     assume each node has already been assigned a globally unique 'name'
#                 
# #             edgeList:
# #                 igraph edge objects
# #                 ref: https://igraph.org/python/doc/igraph.Edge-class.html
#         '''
#         self._graph = graph
#         self._nodes = nodeList
#         add a clusterId to each node # we need to be able to 
#         for each node find the sum of the adj weight and save in dict for quick access
#         
# #         self._edges = edgeList
# #         
# #         self._node2Edges = () # fast look up
# #         self._initNode2Edges()
# #         self._g = ig.Graph()
# #         self._g.add_vertices(nodeList)
# #         self._g.add_edges(edgeList)
#         
#         self._sumOfAdjNodeWeights = {} # key = node['name'], value = weight
#         self._calcSumOfAdjNodeWeights()
#         
#    ############################################################     
#     def _initNode2Edges(self):
#         for e in self._edges:
#                 sadklsdls
#         
#    ############################################################    
#     def _calcSumOfAdjNodeWeights(self):
#         
#    ############################################################    
#     def _findAllEdges(nodeName):    
#         
# #     nodeId = 2 # try and find by name
# #     edges = g.es.select(_source_in = [nodeId])
# #     for e in list(edges):
# #     print(e.tuple)
# #     
# #     edges2 = g.es.select(_target_in = [nodeId])
# #     print("*****")
# #     for e in list(edges2):
# #     print(e.tuple)
# #     
# #     print(type(edges2[0]))
# #     s = [e for e in edges] + [e for e in edges2]
# #     print(len(s))
# #     print(s)        