# Louvain TODO

## input format
- input will be from anndata object
- knn_to_graphModule.py has code convert to igraph
- challenge how are we going to pass this back so that scanpy.pl.louvain can use it
    * [scanpy.tl.louvain](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.louvain.html#scanpy.tl.louvain)
    * from man page
        ```
        Returns
None
By default (copy=False), updates adata with the following fields:

adata.obs['louvain'] (pandas.Series, dtype category)
Array of dim (number of samples) that stores the subgroup id ('0', '1', …) for each cell.

AnnData
        ```