'''
Created on May 16, 2019

@author: andrewdavidson
'''
import logging   
from Edge import Edge
from Node import Node
from Cluster import Cluster


############################################################
class Louvain(object):
    '''
    classdocs
    '''
    logger = logging.getLogger(__name__)

    ############################################################
    def bootStrapInit(self, listOfEdges, listOfWeight):
        '''
        use this constructor to bootstrap from andata.
        
        assigns each cell to its own cluster
        
        ref: knn_to_graphModule get_igraph_from_adjacency()
        
        adjacency = andatq.uns['neighbors']['connectivities']
        sources, targets = adjacency.nonzero()
        listOfWeight = adjacency[sources, targets]
        listOfEdges = list(zip(sources, targets)))
        
        arguments
            listOfEdges: 
                example [(0,2), (2,3)]
                
            listOfWeights:
                example:[ 34.5, 10008.001]
                
        returns a Louvain object
        '''
        self._clusters = []
        
        # TODO: AEDWIP: assert number of nodes == number of cells
        nodesLookup = () # key is nodeId, value is node object
        self._clusterId = 0
        for i in range(len(listOfEdges)):
            # parse input
            edgeTuple = listOfEdges[i]
            weight = listOfWeight[i]
            node1Id, node2Id = edgeTuple
            
            # construct our object
            edge1 = Edge(targetId=node1Id, weight)
            edge2 = Edge(targetId=node2Id, weight)

            self._build(nodesLookup, node1Id, targetEdge=edge2)
            self._build(nodesLookup, node2Id, targetEdge=edge1)

    ############################################################
    def _build(self, nodesLookup, nodeId, targetEdge):
        '''
        todo:
        '''
        if nodeId in nodesLookup:
            n = nodesLookup[nodeId]
        else:
            n = Node(clusterId=self._clusterId, nodeId)
            self._clusterId += 1
            self._nodesLookUp[nodeId] = n
            cluster = Cluster([n])
            self._clusters.append(cluster)
            nodesLookup[nodeId] = n
            
        n.addEdge(targetEdge)

    ############################################################
    def __init__(self, clusters):
        '''
        arguments:
            clusters: a list of cluster objects
        '''
        self._clusters = clusters
        
        self._Q = self._calculateQ()
        
    def _calculateQ(self):
        '''
        calculates modularity
        '''
    
        return Q 