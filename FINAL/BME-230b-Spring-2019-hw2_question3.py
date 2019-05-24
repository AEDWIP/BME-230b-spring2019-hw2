# # BME-230B Spring 2019 HW 3 Question
# James Casaletto, Andrew Davidson, Yuanqing Xue, Jim Zheng
#


from euclid_bbknn import bbknn_graph
import scanpy as sc
import sys
from FStatistic import FStatistic
from bblknn import bblknn_graph
from euclid_knn import KnnG

def main() :
    # check there is exactly one command-line argument provided
    if(len(sys.argv) != 2):
        sys.stderr.write("usage: " + __file__ + " <adata-file-path>\n")
        sys.exit(1)

    # read in adata object from file system
    adata = sc.read(sys.argv[1])

    # build knn graph
    myKnnGraph = KnnG(adata=adata, d_metric='euclidean', n_neighbors=15, method='umap', runPCA=True, nPC=50)

    # run louvain to cluster data
    sc.tl.louvain(myKnnGraph._adata)

    # run umap to project in 2-space
    sc.tl.umap(myKnnGraph._adata)

    # plot the knn graph
    sc.pl.umap(myKnnGraph._adata,  color=['louvain'])

    # ## 3.b. [5 pts]
    # Cluster the integrated dataset using the Louvain method. Re-cluster the data now that you’ve attempted to remove the
    # batch effect. Turn in a UMAP plot showing the integrated dataset and color the cells in the plot by their Louvain
    # cluster assignments.
    #
    # read in ann data file
    anndata = sc.read(sys.argv[1])

    # run our implementation of nearest neighboors and update anndata
    myBBknnGraph = bbknn_graph(adata=anndata, neighbors_within_batch=6, runPCA=True, pcs=50)

    # create louvain clusters
    sc.tl.louvain(myBBknnGraph._adata, flavor='igraph', directed=False, use_weights=True)

    # project data into 2 dimensions
    sc.tl.umap(myBBknnGraph._adata)

    # display graph of louvain clusters
    sc.pl.umap(myBBknnGraph._adata, color=['louvain'])

    # ## 3.c. [10 pts]
    # Quantitatively estimate the degree to which the bb-k-NNG removed the batch
    # effect using the F-statistic described above. Calculate the F statistic using the UMAP
    # solution derived from the original, non-batch balanced 12-k-NNG. Then calculate the F-statistic
    # using the bb-6-NNG to make the UMAP solution. Report both F-statistics. Do you see an
    # improvement in the batch correction using the bb-k-NNG?
    #

    # instantiate and calculate F statistic for non-batch-balanced clusters
    nonbbFstat = FStatistic(myKnnGraph._adata)
    print("non-batch-balanced f stat: " + str(nonbbFstat.calculate_F_statistic()))

    # instantiate and calculate F statistic for batch balanced clusters
    bbFstat = FStatistic(myBBknnGraph._adata)
    print("batch-balanced f stat: " + str(bbFstat.calculate_F_statistic()))

if __name__ == "__main__":
    main()

