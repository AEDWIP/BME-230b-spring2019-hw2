
import igraph
import scanpy as sc

# TEMPLATE
# This is a simplified skeleton following the pseudo code in the assignment.

class Louvain:

    def __init__(self, g, resolution):

        self.g = g
        self.resolution = resolution
    
    def outer_loop(self):
        
        for n_pass in range(self.resolution):
            # Sum of the weights of the links inside each community

            # Sum of the weights of the links incident to nodes in each community
            
            # Loop inner
            self.phase_one()
            
            # Calculate the community set and modularity
            Q = self.modularity()
        
            # Rebuild graph
            # if no community changes then exit outer loop. Possible hint: while loop?
            self.phase_two()
        
        # Update vertices with new communities and edge weights between the two corresponding communities
        
        return communities


    
    def phase_one():

        for vertex_u in community:
            # Find the best community for vertex_u

            if delta_Q > 0:

                # Update the sum of the weights of the links inside each community

                # Update the sum of the weights of the links incident to nodes in each community

                # Update the community information

            # If no vertex moves to a new community then exit inner loop
    
    def modularity():
        Q=0
        for c in community:
            # Compute modularity
        return Q

        
    def phase_two():

        # Update vertices with new community membership

        # Update edges with new community membership. On the first pass this is where we replace the nodes as communities.

        # Sum weights of the links between nodes in the corresponding two communities
        


    def run_louvain(self):

        # Think of structuring communities as a dictionary, where the keys
        # are communities and the values are a list of their nodes
        communities = self.loop_outer()
