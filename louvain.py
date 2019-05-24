'''
Created on May 16, 2019

@author: andrewdavidson
'''
import logging   
from Edge import Edge
from Node import Node
from Cluster import Cluster
from idlelib.idle_test.test_colorizer import source
from scanpy.tools._louvain import louvain


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
    def buildGraph(louvainId, listOfEdges, listOfWeight):
        '''
        use this factory method to bootstrap from anadata.
        
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
        Louvain.logger.info("BEGIN louvainID:{}".format(louvainId))        

        ret = Louvain(louvainId, None)
        ret._leafLouvain = None                
        ret._clusters = dict()  # key is clusterId, value is cluster object
        
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
        
        # TODO: AEDWIP run phase I
        
        Louvain.logger.info("END\n")                
        return ret

   ############################################################
    @staticmethod
    def buildLouvain(louvainId, leafLouvain):
        '''
        use this factory method to buil a new louvain object from the results of a previous run of lovain
        
        implementation starts by constructing a new graph whose nodes are are the clusters
        found in the louvain object argument
        
        returns a louvain object
        '''
        Louvain.logger.info("BEGIN louvainID:{}".format(louvainId))        
        ret = Louvain(louvainId, None)
        ret._leafLouvain = leafLouvain        
        ret._clusters = dict()  # key is clusterId, value is cluster object
        
        # dictionary of all the the nodes in the graph
        # key is nodeId, value is node object
        ret._nodeLookup = {} 
        
        # list of all the edges in the graph
        ret._edges = []
        
        ret._Q = None
        
        ret._phaseII(louvain)
        
        Louvain.logger.info("END\n")        
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
            self._clusters[cluster._clusterId] = cluster
            self._nodeLookup[nodeId] = n
           
        # TODO: AEDWIP: should we be used addEdges ?? 
        n._addEdge(targetEdge) 

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
            if  not nodeI._clusterId == nodeJ._clusterId:
                self.logger.debug("ni:{} cid:{} nj:{} cid:{} not in same cluster"\
                                 .format(nodeI._nodeId, nodeI._clusterId, nodeJ._nodeId, nodeJ._clusterId))
                continue
            else:
                self.logger.debug("ni:{} cid:{} nj:{} cid:{} adding to Q"\
                                 .format(nodeI._nodeId, nodeI._clusterId, nodeJ._nodeId, nodeJ._clusterId))
                
            
            Aij = nodeI.getWeightForEdge(edge._targetId)
            ki = nodeI.getSumAdjWeights()
            kj = nodeJ.getSumAdjWeights()
            i = edge._srcId
            j = edge._targetId
            print('') # AEDWIP:
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
    def changeInModularityIfNodeAdded(self, node, targetCluster):#, isLouvainInit=False
        '''
        calculate change in Q if we add a node to a cluster
        
        formula publish in louvain paper does not work
        '''
        self.logger.debug("BEGIN") 
        
        nodeSet = set()
        if targetCluster._clusterId in node._nodesInClusterDict:
            nodeSet = node._nodesInClusterDict[targetCluster._clusterId]
        else:
#             if isLouvainInit:
#                 # each node is in a separate cluster there are no between edges
#                 return 0
#             else:
            fmt = "target clusterId:{} missing from nodeId:{} _nodesInClusterDict."
            fmt += " should not try  move to clusters node is not connected to "
            self.logger.warn(fmt.format(targetCluster._clusterId, node._nodeId)) 
            return 0
            
        m = self._getM()
        
        ret = 0
        ki = node.getSumAdjWeights()
        for targetNodeId in nodeSet:
            # edges have source and targets
            targetNode = self._nodeLookup[targetNodeId]
            # these edges would no longer in the cluster but between cluster
            kj = targetNode.getSumAdjWeights()
            Aij = node._edgesDict[targetNodeId]._weight
            # multiply by 2 because links are modeled as directed edges
            term = (2 * (Aij - (ki*kj/(2*m))))
            ret += term
            self.logger.debug("ni:{} ti:{} Aij:{} ki:{} kj:{} m:{}"\
                             .format(node._nodeId, targetNodeId, Aij, ki, kj, m))
            
        ret = ret * (1/(2*m))        
        self.logger.debug("END\n") 
        return ret    
    
    ############################################################  
    def changeInModularityIfNodeRemoved(self, node, fromCluster):
        '''
        calculate change in Q if we removed a node from a cluster
        
        formula publish in louvain paper does not work
        '''
        self.logger.debug("BEGIN")  
        
        if not fromCluster._clusterId in node._nodesInClusterDict:
            # there is no requirement that the node is connected to other
            # nodes in the same cluster. I.E. the init step of louvain
            # places each node in a separte cluster
            # we also expect as the algo progresses that some clusters will
            # eventual be emptied
            return 0
        
        nodeSet = node._nodesInClusterDict[fromCluster._clusterId]
        m = self._getM()
        
        ret = 0
        ki = node.getSumAdjWeights()
        for targetNodeId in nodeSet:
            # edges have source and targets
            targetNode = self._nodeLookup[targetNodeId]
            # these edges would no longer in the cluster but between cluster
            kj = targetNode.getSumAdjWeights()
            Aij = node._edgesDict[targetNodeId]._weight
            # multiply by 2 because links are modeled as directed edges
            term = (2 * (Aij - (ki*kj/(2*m))))
            ret += term
            self.logger.debug("ni:{} ti:{} Aij:{} ki:{} kj:{} m:{}"\
                             .format(node._nodeId, targetNodeId, Aij, ki, kj, m))
            
        ret = ret * (1/(2*m))
                
        self.logger.debug("END\n")  
        return ret
        
    ############################################################    
    def _forceAllLazyEval(self):
        for clusterId,c in self._clusters.items():
            c.getSumOfWeights()
            c.getSumOfWeightsInsideCluster(self._nodeLookup)           
    
    ############################################################
    def _getM(self):
        '''
        the m term in the Louvain paper
        "Fast unfolding of communities in large networks"
        
        returns 1/2 the sum of all edges in the graph
        '''
        if not self._m :
            m = 0
            for clusterId, cluster in self._clusters.items():
                m += cluster._getM()
            self._m = m
    
        return self._m
    
    ############################################################
    def getModularity(self): 
        return self.Q   
        
    ############################################################
    def __init__(self, louvainId, clusters=None):
        '''
        TODO"
        calculates modularity 
        
        should only be used by unit test
        
        arguments:
            clusters: a list of cluster objects
            pass None for unit test
        '''
        self._louvainId = louvainId
        self._leafLouvain = None                        
        self._clusters = dict()
        self._Q = None
        self._m = None
        
        # dictionary of all the the nodes in the graph
        # key is nodeId, value is node object
        self._nodeLookup = {} 
        
        # list of all the edges in the graph
        self._edges = []
        
        if not clusters:
            # called from either unit test or buildGraph()
            return   
             
#         if not self._clusters:
#             self.logger.warn("self is not initialized. this is okay if you are running a unit test")
#             return
        
        for c in clusters:
            # TODO: AEDWIP clean up
            self._clusters[c._clusterId] = c        
        
        for clusterId, c in self._clusters.items():
            nodes = c._getNodes()
            for n in nodes:
                if n._nodeId not in self._nodeLookup:
                    self._nodeLookup[n._nodeId] = n
                    self.logger.debug("adding node:{}".format(n._nodeId))
                else:
                    eMsg = "processing cluster:{} node:{} was already in  _nodeLookup".format(c._clusterId, n._nodeId)
                    self.logger.error(eMsg)
                    raise ValueError(eMsg)
                
            self._edges += c._getEdges()

        # TODO: do we need Q? useful for debugging
        self._calculateQ()
                        
    ############################################################  
    def modularityGainIfMove(self, fromCluster, targetCluster, node): # , isLouvainInit=False
        '''
        TODO
        '''     
        self.logger.debug("BEGIN")
        
        changeIfRemoveNode = self.changeInModularityIfNodeRemoved(node, fromCluster)
        changeIfAddNode = self.changeInModularityIfNodeAdded(node, targetCluster)

        ret = changeIfAddNode - changeIfRemoveNode
        self.logger.debug("ret:{} changeIfAddNode:{} loss:{}".format(ret, changeIfAddNode, changeIfRemoveNode, isLouvainInit=False))
        self.logger.debug("END\n")
        return ret
        
    ############################################################ 
    def _phaseI(self, isLouvainInit=False):      
        '''
        TODO
        '''
        self.logger.info("BEGIN")   
          
        self.logger.info("Q:{}".format(self._Q))
        
        # rework initialization
        # make sure cluster is init correctly
        for cId, c in self._clusters.items():
            c.getSumOfWeightsInsideCluster(self._nodeLookup)
            c.getSumOfWeights()
        
        K_CHANGE_IN_Q = 0
        bestMove = (-1, -1, -1, -1) # (changeInQ, node, fromCluster, toCluster
        isImproving = True
        while isImproving:
            self.logger.info("BEGIN EPOCH")
            isImproving = False
            
            # links in graph are modeled as a pair of directed edges
            # keep track of moves so we do not cycle
            # i.e. move na from cluster 1 to cluster 2, then back again
            # use a set of tuples. each element is of form (srcNodeId, targetNodeId)
            trackMoves = set()
                            
            for nodeId, node in self._nodeLookup.items():
                # keep track of which clusters we have tested
                testedClusters = set()
                
                # Q will only improve if we are moving into a cluster that has a node 
                # we are connected to
                fromCluster = self._clusters[node._clusterId]
                for candidateClusterId, candidateNodeSet in node._nodesInClusterDict.items():
                    if node._clusterId == candidateClusterId:
                        continue
                    
                    #targetCluster = self._clusters[candidateClusterId]
                    for candidateNodeId in candidateNodeSet:
                        candidateNode = self._nodeLookup[candidateNodeId]
                        
                        targetClusterId = candidateNode._clusterId
                        if fromCluster._clusterId == targetClusterId:
                            continue
                        
                        targetCluster = self._clusters[targetClusterId]
                        if targetCluster in testedClusters :
                            continue
                        testedClusters.add(targetCluster)
                        
                        # track moves to prevent cycles
                        possibleMove = (nodeId, targetClusterId)
                        if possibleMove in trackMoves :
                            continue
                        trackMoves.add(possibleMove)              
                        
                        predictedChange = self.modularityGainIfMove(fromCluster, targetCluster, node)
                        if predictedChange > bestMove[K_CHANGE_IN_Q] and predictedChange > 0:
                            bestMove = (predictedChange, node, fromCluster, targetCluster)
                        
                if bestMove[0] > 0:
                    isImproving = True
                    change, node, fromC, toC = bestMove
                    bestMove = (-1, -1, -1, -1) # (changeInQ, node, fromCluster, toCluster

                    self._Q += change
                    print('')
                    self.logger.info("Q:{}, change:{} nodeId:{} fromClusterId:{} toClusterId:{}"\
                                     .format(self._Q, change, node._nodeId, fromC._clusterId, toC._clusterId))
                    fromCluster.moveNode(targetCluster, node, self._nodeLookup, isLouvainInit)
                    
            # TODO: prune empty clusters
                    
            print('')
            for cid,c in self._clusters.items():
                self.logger.info(c)
            self.logger.info("END of epoch\n")
                
            
        self.logger.info("Q:{}".format(self._Q))        
        self.logger.info("END\n")     
        
        
    ############################################################ 
    def _phaseIICreateNewEdges(self, isLouvainInit=False):      
        '''
        TODO
        creates edges from leafLouvain clusters
        
        does not calculate weights use betweenEdgeWeightsDict to adjust weights
        
        returns (nodeEdgesDict, betweenEdgeWeightsDict)
            nodeEdgesDict 
                key is new nodeId, value set of edges
            
            betweenEdgeWeightsDict 
                key is (srcId:targetId) value is list of weight of leaf edges
                be careful not to double weight. keys will be double entry (a,b) and (b,a)
        '''
        self.logger.info("BEGIN")   
        
        # key is new nodeId == leaf clusterId: set of edges
        # use set to prevent duplicate edges
        nodeEdgesDict = dict() 
        UNKNONWN_WEIGHT = None # place holder. we need to calculate weights later
        
        # key is tuple (srcId:targetId) value is list of weight of leaf edges
        # use a list to make debugging easier
        # most of our unit test use edge weight = 1. anadata will be distance value
        # be careful not to double weight. keys will be double entry (a,b) and (b,a)
        betweenEdgeWeightsDict = dict() 
        
        for leafNodeId, leafNode in self._leafLouvain._nodeLookup.items():
            leafNodeClusterId = leafNode._clusterId
            
            # to make debugging easier our new node ids will be the leaf cluster ids            
            newNodeId = leafNodeClusterId
            
            for adjLeafNodeClusteId, adjLeafNodeSet in leafNode._nodesInClusterDict.items():
                if adjLeafNodeClusteId == leafNodeClusterId:
                    # edges inside leaf clusters do not create edges between cluster
                    # in phase II.  
                    continue
               
                # create edges between clusters
                if adjLeafNodeSet:
                    # https://www.pythoncentral.io/how-to-check-if-a-list-tuple-or-dictionary-is-empty-in-python/
                    # weird python syntax to test if list, set, dict is empty or not
                    
                    if not newNodeId in nodeEdgesDict:
                        nodeEdgesDict[newNodeId] = set()
                        
                    for adjNodeId in adjLeafNodeSet:
                        adjNode = self._leafLouvain._nodeLookup[adjNodeId]
                        newTargetNodeId = adjNode._clusterId
                        e = Edge(weight=UNKNONWN_WEIGHT, srcId=newNodeId, targetId=newTargetNodeId)
                        nodeEdgesDict[newNodeId].add(e) 
                        
                        if not leafNodeClusterId in betweenEdgeWeightsDict:
                            betweenEdgeWeightsDict[newNodeId] = []
                        
                        w = adjNode._edgesDict[leafNodeId]._weight
                        key = (newNodeId,newTargetNodeId)
                        if not key in betweenEdgeWeightsDict:
                            betweenEdgeWeightsDict[key] = []
                        betweenEdgeWeightsDict[key].append(w)
                        
        self.logger.info("END\n") 
        return nodeEdgesDict, betweenEdgeWeightsDict, 
        

    ############################################################     
    def _fixThisBugUseTrueOOEncapsliations(self, nodeEdgesDict, betweenEdgeWeightsDict):
        '''
        do not be lazy! it leads to bugs
        
        coded in haste. assume oh all the code is same package
        no need to implement accessor fuctions. This lead to 
        lots of issue.
        
        we need to construct the full object graph before calculating
        any of the values
        '''
        # calculate edge weights
        for nodeId in nodeEdgesDict.keys():
            edges = nodeEdgesDict[nodeId]
            for e in edges:
                key = (nodeId,e._targetId)
                listOfWeights = betweenEdgeWeightsDict[key]
                e._weight = sum(listOfWeights)
                
        print()
        for k,v in nodeEdgesDict.items():
            self.logger.info("nodeEdgesDict nodeId:{} edges:{}".format(k,v))
            
        print()
        for k,v in betweenEdgeWeightsDict.items():
            self.logger.info("betweenEdgeWeightsDict key:{} listOfWeights:{}".format(k,v))        
                            
        # create nodes and clusters
        for newNodeId in nodeEdgesDict.keys():
            # create new node
            newClusterId = newNodeId
            newNode = Node(newClusterId, newNodeId) 
            self._nodeLookup[newNodeId] = newNode
            # create new cluster
            newCluster = Cluster(newClusterId, [newNode])
            self._clusters[newClusterId] = newCluster            
            
        # add edges to nodes       
        for newNodeId, edgeSet in nodeEdgesDict.items():
            #edgeList = nodeEdgesDict[newNodeId]
            #newNode.addEdges(edgeList,  self._nodeLookup) 
            newNode = self._nodeLookup[newNodeId]
            for e in edgeSet: 
                newNode._addEdge(e)
                
        
        # init node caches
        for nId in self._nodeLookup.keys():
            node = self._nodeLookup[nId]
            # because we used _addEdge() instead of addEdges()
            # we need to make sure cache is set up
            node._initKiinCache(self._nodeLookup)      
            
        # force nodes to calc cached values
        for nodeId in self._nodeLookup.keys():
            node = self._nodeLookup[nodeId]
            node.getSumAdjWeights()
            node.getSumOfWeightsInsideCluster(nodeId, self._nodeLookup)
            
        # force clusters to calc cached values
        for clusterId in self._clusters.keys():
            # run lazy eval
            cluster = self._clusters[clusterId]
            cluster.getSumOfWeights()
            cluster.getSumOfWeightsInsideCluster(self._nodeLookup) 
                               
        self.logger.info("END\n")     
    
    def _phaseII(self, isLouvainInit=False):      
        '''
        TODO creates graph from leafLouvain
        '''
        self.logger.info("BEGIN") 
        nodeEdgesDict, betweenEdgeWeightsDict = self._phaseIICreateNewEdges()    
        self._fixThisBugUseTrueOOEncapsliations(nodeEdgesDict, betweenEdgeWeightsDict) 
        
        
    ############################################################                
    def __repr__(self):
        self._forceAllLazyEval()
        ret = "\n\tid:{}".format(self._louvainId)
        ret += "\tQ:{}".format(self._Q)
        ret += "\tnumber of Nodes:{}\n".format(len(self._nodeLookup.keys()))
        ret += "\tnumber of edges:{}\n".format(len(self._edges))
        ret += "\tnumber of clusters:{}\n".format(len(self._clusters.keys()))
        for c in self._clusters:
            ret += "\t{}\n".format(c)
                        
        return ret        
