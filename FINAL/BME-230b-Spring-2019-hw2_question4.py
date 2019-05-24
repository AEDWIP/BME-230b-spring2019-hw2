# # BME-230B Spring 2019 HW 4 Question
# James Casaletto, Andrew Davidson, Yuanqing Xue, Jim Zheng
#

# ## 4.b. [10 pts]
# Turn in a bar plot of the Adjusted Rand Index (ARI) for Louvain clusters obtained from 10 independently subsampled
# bb-3-6-NNGs compared to the Louvain clusters obtained on the original bb-6-NNGs. Also report the average and standard
# deviations of the ARI. Based on these results, would you conclude these clusters are robust? Justify your answer.
# Hint: check if your ARI is significantly better than chance.
#

from euclid_bbknn import bbknn_graph
import matplotlib.pyplot as plt
import sys
import scanpy as sc
from ARIStatistic import ARIStatistic
import numpy as np


def main() :
    # check there is exactly two command-line arguments provided
    if(len(sys.argv) != 3):
        sys.stderr.write("usage: " + __file__ + " <adata-file-path> <number-of-samples>\n")
        sys.exit(1)

    # read in adata object from file system
    adata_bbknn = sc.read(sys.argv[1])

    # run bbknn on adata
    bbknn = bbknn_graph(adata_bbknn, neighbors_within_batch=6, runPCA=True, pcs=50)

    # run louvain to cluster data
    sc.tl.louvain(bbknn._adata)

    # read in AnnData object (again)
    adata_blknn = sc.read(sys.argv[1])

    # instantiate ARIStatistic object
    myARIStat = ARIStatistic(adata=adata_blknn, k_per_batch=6, l=3, n_components=50, n_samples=int(sys.argv[2]),
                             bbknn_louvain=bbknn._adata.obs['louvain'])

    # plot ARIs
    plt.bar(range(1, len(myARIStat._results) + 1), myARIStat._results, color='black')

    # print statistics
    print('Average:', np.mean(myARIStat._results), 'SD:', np.std(myARIStat._results))

if __name__ == "__main__":
    main()





