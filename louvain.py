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
    public functions:
    
        bootStrapInit(self, listOfEdges, listOfWeight)
        getModularity()
        __init__(self, clusters)
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
        
        # dictionary of all the the nodes in the graph
        # key is nodeId, value is node object
        self._nodeLookup = {} 
        
        # list of all the edges in the garpy
        self._edges = []
        
        self._Q = None
        
        # TODO: AEDWIP: assert number of nodes == number of cells
        self._clusterId = 0
        for i in range(len(listOfEdges)):
            # parse input
            edgeTuple = listOfEdges[i]
            weight = listOfWeight[i]
            node1Id, node2Id = edgeTuple
            
            # construct our object
            edge1 = Edge(weight, srcId=node1Id, targetId=node2Id)
            edge2 = Edge(weight, srcId=node2Id, targetId=node1Id)
            
            # edge1 and edge 2 are identical 
            self._edges.append(edge1)
            # TODO:AEDWIP double check if we need to either remove 1/2 or add edge2

            self._build(node1Id, targetEdge=edge1)
            self._build(node2Id, targetEdge=edge2)
            
        self._calculateQ()

    ############################################################
    def _build(self, nodeId, targetEdge):
        '''
        todo:
        '''
        if nodeId in self._nodeLookup:
            n = self._nodeLookup[nodeId]
        else:
            n = Node(self._clusterId, nodeId)
            self._nodeLookup[nodeId] = n
            cluster = Cluster(self._clusterId, [n])
            self._clusterId += 1            
            self._clusters.append(cluster)
            self._nodeLookup[nodeId] = n
            
        n.addEdge(targetEdge)

    ############################################################
    def __init__(self, clusters):
        '''
        arguments:
            clusters: a list of cluster objects
        '''
        eMsg = "AEDWIP NOT IMPLEMENTED YET!"
        self.logger.error(eMsg)
        
    
    ############################################################
    def _calculateM(self):
        '''
        the m term in the Louvain paper
        "Fast unfolding of communities in large networks"
        
        returns 1/2 the sum of all edges in the graph
        '''
        
        m = 0
        for cluster in self._clusters:
            m += cluster._calculateM()
    
        return m 
    
    ############################################################
    def _calculateQ(self):
        '''
        calculates modularity for graph
        
        from igraph? looks wrong?
        
        Q=1/(2m) * sum( Aij-ki*kj / (2m) delta(ci,cj),i,j). 
        m is the number of edges, 
        Aij is the element of the A adjacency matrix in row i and column j, 
        ki is the degree of node i, 
        kj is the degree of node j, 
        Ci and cj are the types of the two vertices (i and j). 
        delta(x,y) is one iff x=y, 0 otherwise.
        '''
        self.logger.info("BEGIN")
        # note implemenation already multiplied by 1/2
        # TODO should we move the mulply out? it should technically be faster
        # does not effect Big O
        m = self._calculateM()
        self.logger.info("m:{}".format(m))
        
        modularitySumTerm = 0
        for edge in self._edges:
            self.logger.info("\n{}".format(edge))
            
            nodeI = self._nodeLookup[edge._srcId]
            nodeJ = self._nodeLookup[edge._targetId]
            self.logger.info("nodeI:{}".format(nodeI))
            self.logger.info("nodeJ:{}".format(nodeJ))
            
            # calculate the sigma term
            if not (nodeI._clusterId == nodeJ._clusterId):
                self.logger.info("not in same cluster")
                continue
            
            Aij = nodeI.getWeightForEdge(edge._targetId)
            ki = nodeI.getSumAdjWeight()
            kj = nodeI.getSumAdjWeight()
            term = Aij - ki*kj / 2*m
            self.logger("(Aij:{} - ki:{}*kj:{}/2m:{}) == {}".format(Aij, ki, kj, m, term))
            modularitySumTerm += term 
        

        self._Q = modularitySumTerm/m
        self.logger.info("END\n")

    
    ############################################################
    def getModularity(self): 
        return self.Q   