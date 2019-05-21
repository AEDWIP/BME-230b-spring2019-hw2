
import igraph
import scanpy as sc


class Louvain:
    """
    Fill in class description!
    """

    def __init__(self, g, adata):

        # iGraph object
        self.g = g
        self.adata = adata
    
    def assign_community_id(self):
        """
        Assign each vertex a commmunity id
        """
        # Note: this method does not have to return anything
        # Hint: read test_graph()
        return None

    def create_node_to_community_tracker(self):
        """
        Create a data structure to keep track of nodes and their community id's
        """
        return community
    
    def create_adjacency_structure(self):
        """
        Create adjacency table structure
        """
        return adjacency
    
    def compute_m(self):
        """
        Sum of the weights of all the links in the network
        """
        return m
    
    def phase_one(m, community, adjacency):
        """
        Fill in doc string!
        """
        # Write loop inner and sections 3.1-5.4
        # You will have to rearrange code!
        # Consider setting up a test for phase one

        # Section 3.1
        for node_i in community:

            # Find all node_i's neighbor(s) j

            # Section 3.2
            # Remove i from its community and place it in community j

            ### Section 4: Find the best community for node_i ###
            # Note: you must determine what parameters to pass
            # to each method in section 4.1-4.4

            e_in = self.compute_sigma_in() # Section 4.1
            e_tot = self.compute_sigma_total() # Section 4.2
            ki = self.compute_ki() # Section 4.3
            ki_in = self.compute_ki_in() # Section 4.4

            # Compute delta Q
            delta_q = self.compute_delta_q(e_in, ki_in, e_tot, ki, m) # Section 4.5

            ### Section 5 ###
            if delta_q > 0: # Section 5.1
                # Section 5.2: Edge case alert!
                # Consider how you might make sure
                # vertex u finds the best community
                # if another vertex j's community
                # also has a positive modularity.

                # Section 5.3
                # Update node_i's community id
        
        # Section 5.4
        # If no vertex moves to a new community then exit inner loop.
        # Note that this section of code will need to be written outside 
        # of section 3.1

        # Hint: What data structure can you create that can help you keep
        # track of the sum of the weights of the links incident to node_i? 
        # In a second pass these nodes will be communities!
        # I've named this weight_sum_incident_to_node but feel free to change

        return weight_sum_incident_to_node, community
    
    def compute_sigma_in(self):
        """
        Fill in doc string!
        """
        return e_in
    
    def compute_sigma_total(self):
        """
        Fill in doc string!
        """
        return e_tot
    
    def compute_ki(self):
        """
        Fill in doc string!
        """
        return ki
    
    def compute_ki_in(self):
        """
        Fill in doc string!
        """
        return ki_in
    
    def compute_delta_q(self, e_in, ki_in, e_tot, ki, m):
        """
        Fill in doc string!
        """

        delta_q = ((e_in + ki_in)/ (2*m)  - ((e_tot + ki)/ (2*m)) ** 2) - (((e_in)/(2*m)) - ((e_tot)/(2*m)) **2 - ((ki)/(2*m))**2)

        return delta_q
    
    def compute_modularity(self):
        """
        Fill in doc string!
        """
        return modularity

    def phase_two(self):
        """
        Fill in doc string!
        """
        # Section 7.1: Reassign nodes with new community ids
        # Section 7.2: Update edge weights with new community id's

        # Hint: Reference images in section 7.2.1 and 7.2.2 may be helpful to you
        return updated_community, updated_adjacency
    
    def update_adata(self):
        """
        Fill in doc string!
        """
        return updated_adata
    
    def run(self):
        # Note sections and subsections are referenced to notes uploaded on canvas
        
        ### Section 1: Initialize data structures and community id ###
        self.assign_community_id() # Section 1.1
        community = self.create_node_to_community_tracker() # Section 1.2
        adjacency = self.create_adjacency_structure() # Section 1.3
        m = self.compute_m() # Section 1.4
        
        ### Section 2: Loop Outer  ###
        # Fill in section. You will have to rearrange code!
        # Check yourself before you wreck yourself: 
        # 1) Do you need to implement the entire algorithm at once? 
        # 2) How can you split this up into sections and test phase one then two?
        
        ### Section 3: Phase One and Loop Inner ###
        weight_sum_incident_to_node, community = self.phase_one(m, community, adjacency)

        # Section 6: Compute community set and modularity
        # Note: you must determine what parameters to pass to this method
        modularity = self.compute_modularity() 

        ### Section 7: Rebuild Graph ###
        # Note: you must determine what parameters to pass to this method
        updated_community, updated_adjacency = self.phase_two()

        ### Section 8 ###
        # If no community changes then exit outer loop and return updated graph
        # and it's modularity score

        # Update anndata object in the observation data frame with new cluster id's
        # Note: you must determine what parameters to pass to this method
        updated_adata = self.update_adata()

        return updated_adata, modularity

def test_graph():
    """
    Sample graph as illustrated in notes for testing your algorithm!
    """
    
    # Create sample graph
    sample_graph = [(0,1), (0,2), (2,3), (3,4), (4,2), (2,5), (5,0), (6,3), (5,6)]
    g = ig.Graph(sample_graph)
    
    # Assign each vertex a community id
    # Note: You will have to implement a similar function in the Louvain class 
    # Answer bite: g.vs['community'] = [x for x in range(len(sample_graph)]
    
    # In the test graph I manually assign a weight to each edge.
    # However, you will be using the umap connectivities
    weight = 1
    for edge in g.es:
        edge['weight'] = weight
        weight +=1

    return g

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
            logg.warn('The constructed graph has only {} nodes. '
                      'Your adjacency matrix contained redundant nodes.'
                      .format(g.vcount()))
        return g

def main():

    # Remember you must use your KNN 
    # adata = sc.read('PBMC.merged.h5ad')
    # adjacency = adata.uns['neighbors']['connectivities']
    # g = get_igraph_from_adjacency(adjacency, directed=False)
    
    adata, modularity = Louvain.run(g, adata)

if __name__ == '__main__':
    main()