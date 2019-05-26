# # BME-230B Spring 2019 HW 2 Question 1
# James Casaletto, Andrew Davidson, Yuanqing Xue, Jim Zheng
# 
# Question 1.a, 1.b see [euclid_knn.py](euclid_knn.py)
# 
# ref: 
# - [ Single-Cell Analysis in Python](https://scanpy.readthedocs.io/en/stable/api/index.html#tools-tl)
# - [data exploration](exploreData.ipynb)
# - [scanpy.tl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.umap.html)
# - <span style="color:red">scanpy.api.pl no longer exists</span>
# - [scanpy.pl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.pl.umap.html#scanpy.pl.umap)
# 

from euclid_knn import KnnG
import matplotlib.pyplot as plt
import numpy as np
import scanpy.api as sc
import scanpy
import sys


def main() :
    # check there is exactly one command-line argument provided
    if(len(sys.argv) != 2):
        sys.stderr.write("usage: " + __file__ + " <adata-file-path>\n")
        sys.exit(1)

    anndata = sc.read(sys.argv[1])

    # run PCA, compute pairwise distance, k-nearest-neighbors \n# trick to replacing scanpy implementation with our own
    # is to \n# update the anndata object with our intermediate values\n# \n# we replace the scanpy version of PCA by
    # updating\n# .adata.obsm['X_pca'] = our PCA() output\n#\n# we replace the scanpy version of k-nearest-neighbors by
    # updating\n# self.adata.uns['neighbors']['connectivities'] = our knn() output\n
    # self.adata.uns['neighbors']['distances'] = out knn()
    # output\n\nknng = KnnG(anndata, n_neighbors=12, runPCA=True, nPC=50)")

    # umap() reduce results to 2 deminsions so that we can plot the data\nsc.tl.umap(anndata)')


    # ### 1.c. [5 pts] Turn in a UMAP plot of your 12-NN graph calculated from the combined chemistry PBMC dataset
    # colored by batch (the chemistry used)
    scanpy.pl.umap(anndata, color=['Method'])

    # ### 1.d. [5 pts] Turn in another UMAP plot of your 12-NN graph calculated from the combined chemistry PBMC dataset
    # but colored by cell type
    scanpy.pl.umap(anndata, color=['Cell type'])

if __name__ == "__main__":
    main()
