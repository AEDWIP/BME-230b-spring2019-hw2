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
    
        @staticmethod
        buildGraph(self, listOfEdges, listOfWeight)
        
        __init__(self, clusters)

        getModularity()
    '''
    logger = logging.getLogger(__name__)

    ############################################################
    @staticmethod
    def buildGraph(listOfEdges, listOfWeight):
        '''
        use this constructor to bootstrap from andata.
        
        assigns each cell to its own cluster and calculates modularity
        
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
        ret = Louvain(None)
        ret._clusters = []
        
        # dictionary of all the the nodes in the graph
        # key is nodeId, value is node object
        ret._nodeLookup = {} 
        
        # list of all the edges in the graph
        ret._edges = []
        
        ret._Q = None
        
        # TODO: AEDWIP: assert number of nodes == number of cells
        
        ret._clusterId = 0
        tmpDirectedEdgeLookup = set() # elements are (node1Id,node2Id)
        for i in range(len(listOfEdges)):
            # parse input
            edgeTuple = listOfEdges[i]
            weight = listOfWeight[i]
            node1Id, node2Id = edgeTuple
            
            #
            # construct our object
            # our design is node centric. To make debug, and unit test easier
            # we want our nodes to have easy access to all their adj edges
            # so if we have and edge from b to c, we want to insure there is an edge 
            # from c to b. If the listOfEdges was generated from a pair wise distance
            # matrix it would automatically contain both edge. How eve we are working
            # with knn edges. given b is a nearest neighbor of c, does not imply c is
            # nearest neighbor of b
            #
            # if b and c are symmetric, ie. b and c are both nearest neighbors of each other
            # we need to be careful we do not add the same edge to the graph multiple times
            #
            
            if not (node1Id, node2Id) in tmpDirectedEdgeLookup:
                tmpDirectedEdgeLookup.add((node1Id, node2Id))
                edge1 = Edge(weight, srcId=node1Id, targetId=node2Id)
                ret._edges.append(edge1)
                ret._build(node1Id, targetEdge=edge1)

            if not (node2Id, node1Id) in tmpDirectedEdgeLookup:
                tmpDirectedEdgeLookup.add((node2Id, node1Id))  
                edge2 = Edge(weight, srcId=node2Id, targetId=node1Id)
                ret._edges.append(edge2)
                ret._build(node2Id, targetEdge=edge2)
                
        for nodeId, node in ret._nodeLookup.items():
            node._initKiinCache(graphNodesLookup=ret._nodeLookup)

        ret._calculateQ()
        return ret

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
            
        n._addEdge(targetEdge) 

    ############################################################
    def __init__(self, clusters=None):
        '''
        TODO"
        calculates modularity 
        
        should only be used by unit test
        
        arguments:
            clusters: a list of cluster objects
            pass None for unit test
        '''
        
        if not clusters:
            # called from either unit test or buildGraph()
            return
            
        self._clusters = clusters
        self._Q = None
        self._m = None
        
        # dictionary of all the the nodes in the graph
        # key is nodeId, value is node object
        self._nodeLookup = {} 
        
        # list of all the edges in the graph
        self._edges = []
        
        if not self._clusters:
            self.logger.warn("self is not initialized. this is okay if you are running a unit test")
            return
        
        for c in self._clusters:
            nodes = c._getNodes()
            for n in nodes:
                if n._nodeId not in self._nodeLookup:
                    self._nodeLookup[n._nodeId] = n
                    self.logger.info("adding node:{}".format(n._nodeId))
                else:
                    eMsg = "processing cluster:{} node:{} was already in  _nodeLookup".format(c._clusterId, n._nodeId)
                    self.logger.error(eMsg)
                    raise ValueError(eMsg)
                
            self._edges += c._getEdges()

        
        # TODO: do we need Q? useful for debugging
        self._calculateQ()
        
    
    ############################################################
    def _getM(self):
        '''
        the m term in the Louvain paper
        "Fast unfolding of communities in large networks"
        
        returns 1/2 the sum of all edges in the graph
        '''
        
        if not self._m :
            m = 0
            for cluster in self._clusters:
                m += cluster._getM()
            self._m = m
    
        return self._m
    
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
        self.logger.debug("BEGIN")
        # note implemenation already multiplied by 1/2
        # TODO should we move the mulply out? it should technically be faster
        # does not effect Big O
        m = self._getM()
        self.logger.debug("m:{}".format(m))
        
        modularitySumTerm = 0
        for edge in self._edges:
            self.logger.debug("\n{}".format(edge))
            
            nodeI = self._nodeLookup[edge._srcId]
            nodeJ = self._nodeLookup[edge._targetId]
            self.logger.debug("nodeI:{}".format(nodeI))
            self.logger.debug("nodeJ:{}".format(nodeJ))
            
            # calculate the sigma term
            if not (nodeI._clusterId == nodeJ._clusterId):
                self.logger.debug("not in same cluster")
                continue
            
            Aij = nodeI.getWeightForEdge(edge._targetId)
            ki = nodeI.getSumAdjWeights()
            kj = nodeJ.getSumAdjWeights()
            i = edge._srcId
            j = edge._targetId
            self.logger.debug("i:{} j:{} A{}{}:{} k{}:{} k{}:{} 2*m:{}".format(i,j, i, j, Aij, i, ki, j, kj, 2*m))
            self.logger.debug(" ki*kj / 2*m == {}".format( (ki*kj) / (2*m)))
            term = Aij - (ki*kj) / (2*m) 
            self.logger.debug("(Aij:{} - ki:{}*kj:{}/2m:{}) == {}".format(Aij, ki, kj, m, term))
            modularitySumTerm += term 
        

        self.logger.debug("m:{} modularitySumTerm:{}".format(m, modularitySumTerm))
        self._Q = modularitySumTerm/(2*m)
        
        if not (-1.0 <= self._Q <= 1.0):
            eMsg = "invalid Q=={} must be -1.0 <= Q <= 1.0".format(self._Q)
            self.logger.error(eMsg)
            raise ValueError(eMsg)
        
        self.logger.debug("END\n")
    
    ############################################################
    def getModularity(self): 
        return self.Q   
    
    ############################################################    
    def _forceAllLazyEval(self):
        for c in self._clusters:
            c.getSumOfWeights()
            c.getSumOfWeightsInsideCluster(self._nodeLookup)   
        
    ############################################################                
    def __repr__(self):
        self._forceAllLazyEval()
        ret = "\n\tQ:{}".format(self._Q)
        ret += "\tnumber of Nodes:{}\n".format(len(self._nodeLookup.keys()))
        ret += "\tnumber of edges:{}\n".format(len(self._edges))
        ret += "\tnumber of clusters:{}\n".format(len(self._clusters))
        for c in self._clusters:
            ret += "\t{}\n".format(c)
            
            
        return ret
    

    ############################################################  
    def modularityGainIfMove(self, targetCluster, node, graphNodesLookup): 
        '''
        TODO
        '''     
        self.logger.info("BEGIN")
        
        # calculate change in Q cause by removing node
        sigmaIn = targetCluster.getSumOfWeightsInsideCluster(self._nodeLookup)
        kiin = node.getSumOfWeightsInsideCluster(targetCluster._clusterId, self._nodeLookup)
        m = self._getM()
        sigmaTot = targetCluster.getSumOfWeights()
        ki = node.self.getSumAdjWeights()
        
        loss = ((sigmaIn - kiin)/(2*m)) + ((sigmaTot - ki)/(2*m))**2

        
        # calculate change Q caused by adding node
        gain = ((sigmaIn + kiin)/(2*m)) - ((sigmaTot + ki)/(2*m))**2

        self.logger.info("END\n")
        return gain + loss
        
