import igraph as ig
import numpy as np
import networkx as nx
import scanpy as sc
import collections

# JC

def get_igraph_from_adjacency(adjacency, directed=None):
        """Get igraph graph from adjacency matrix."""
        sources, targets = adjacency.nonzero()
        weights = adjacency[sources, targets]
        if isinstance(weights, np.matrix):
            weights = weights.A1
        g = ig.Graph(directed=directed)
        g.add_vertices(adjacency.shape[0])  # this adds adjacency.shap[0] vertices
        g.add_edges(list(zip(sources, targets)))
        try:
            g.es['weight'] = weights
        except:
            pass
        if g.vcount() != adjacency.shape[0]:
            log.warn('The constructed graph has only {} nodes. '
                      'Your adjacency matrix contained redundant nodes.'
                      .format(g.vcount()))
        return g

def convertIGraphToNxGraph(myig):
    node_names = myig.vs["community"]
    edge_list = myig.get_edgelist()
    weight_list = myig.es["weight"]
    node_dict = collections.defaultdict(str)

    for idx, node in enumerate(myig.vs):
        node_dict[node.index] = node_names[idx]

    convert_list = []
    for idx in range(len(edge_list)):
        edge = edge_list[idx]
        new_edge = (node_dict[edge[0]], node_dict[edge[1]], weight_list[idx])
        convert_list.append(new_edge)

    convert_graph = nx.Graph()
    convert_graph.add_weighted_edges_from(convert_list)
    return convert_graph

pbmc = sc.read('/Users/jcasaletto/PycharmProjects/BME230B/HW2/BME-230b-spring2019-hw2/PBMC.merged.h5ad')
adjacency = pbmc.uns['neighbors']['connectivities']

# get igraph from adjacency matrix
g = get_igraph_from_adjacency(adjacency, directed=False)

# Assign community labels to each node
g.vs['community'] = [x for x in range(pbmc.shape[0])]

# xNetwork
xg = convertIGraphToNxGraph(g)
