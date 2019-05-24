# # BME-230B Spring 2019 HW 2 Question
# James Casaletto, Andrew Davidson, Yuanqing Xue, Jim Zheng
# 
# - ref
#  * [scanpy.tl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.umap.html)
#  * [scanpy.api.pp.neighbors](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.api.pp.neighbors.html?highlight=neighbors)
#  * [scanpy.pl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.pl.umap.html#scanpy.pl.umap)
#  * [scanpy.tl.louvain](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.louvain.html#scanpy.tl.louvain)
#  * [GSEAPY: Gene Set Enrichment Analysis in Python. pypi.org](https://pypi.org/project/gseapy/)
#  * [GSEAPY: Gene Set Enrichment Analysis in Python gseapy.readthedocs.io](https://gseapy.readthedocs.io/en/latest/introduction.html)
#  * [anndata](https://anndata.readthedocs.io/en/latest/anndata.AnnData.html)
#         + "uns" stands for unstructured data
#         + "obs" are panda data frame observations 
#         + "obsm key-indexed multi-dimensional observations
#  * [Hypergeometric_distribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution)
#  * [Hypergeometric Tests
# for Gene Lists](http://users.unimi.it/marray/2007/material/day4/Lecture7.pdf)



from euclid_knn import KnnG
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy.api as sc
import scanpy
import scipy.special
import scipy.stats as stats
import sys


def main() :
    # check there is exactly one command-line argument provided
    if(len(sys.argv) != 2):
        sys.stderr.write("usage: " + __file__ + " <adata-file-path>\n")
        sys.exit(1)

    # read in AnnData from file
    anndata = sc.read(sys.argv[1])

    # ## 2.b. [5 pts]
    # Turn in a UMAP plot of the combined dataset as you did in question #1, but
    # this time, color the cells by their Louvain cluster assignments determined for each cell
    # within each batch as a different color in each plot.

    # run our implementation of nearest neighboors and update anndata
    KnnG(anndata, n_neighbors=12, runPCA=True, nPC=50)

    # running Scanpy's version of Louvian
    scanpy.tl.louvain(anndata, flavor='igraph', directed=False, use_weights=True)

    # umap the coordinates to 2-space
    scanpy.tl.umap(anndata)

    # display clustering by cell type, louvain clustering, and batch
    sc.pl.umap(anndata,  color=['Cell type', 'batch','louvain'])

    # ## 2.c. [5 pts]
    # Turn in a table that lists each cluster and its best-matching cell type
    # annotation. The table should contain the cluster number and its best matching cell-type
    # annotation based on the hypergeometric analysis.

    # calculate counts of cells and cell types per cluster
    cellCountsByClusterId, cellTypesInClusters = hw2q2.createCountsDict(anndata)

    # find best annotation for cells in each cluster using hypergeometric p-value
    annotations = hw2q2.bestAnnotation(anndata, cellTypesInClusters, cellCountsByClusterId)

    # sort output by clusterId (which is ascii string, so we need to first create an integer version
    # sort on the integer column, then drop the integer column before displaying data frame
    annotations['clusterIdNumber'] = [int(i) for i in annotations['clusterId']]
    annotations.sort_values(by=['clusterIdNumber'], ascending=[True], inplace=True)
    annotations.drop(columns=['clusterIdNumber'], axis=1, inplace=True)
    print(annotations)

    # ## 2.d. [5 pts]
    # Turn in a list of top 5 pathways for each cluster in each dataset. You should
    # use the gene expression signature of each cluster to find an associated pathway. A gene
    # signature for a cluster represents the gene expression levels for a characteristic cell that is a
    # member of the cluster. Use the centroid 𝞵 i of the i th cluster as the signature. Compute the
    # centroids for each cluster in each dataset. You will next derive a gene-signature based
    # annotation for each cluster using these centroids. Use a list of Gene Ontology Biological
    # Process categories (provided in the Resources section at the top of this homework) and your
    # signatures to perform an all-against-all Gene Set Enrichment Analysis (GSEA). Turn in a table
    # that lists the top 5 pathways for each cluster

    # 1. create a data frame we can use to select the cells in a given cluster
    #     a. we need the index value we can use to to get the cell's expression
    #     values from the numpy array anndata.obsm['X_pca'].shape
    #testGetCellsIdxForCluster(anndata)
    #testGetGeneExpressionSignatureForCluster(anndata)

    # find centroids of clusters
    clusterSigs = hw2q2.calculateClusterSignatures(anndata)

    # use gseapy to provide top 5 pathways for each cluster
    pathways = hw2q2.rankPathWays(anndata, clusterSigs, topN=5)

    # print pathways
    print(pathways)

class hw2q2():
    @staticmethod
    def createCountsDict(anndata):
        '''
        assumes the data has been grouped by cluster id and cell type.

        return:
            cellCountsByClusterId:
                a dictionary of form d['cellType']['cellType'] == count

            cellTypesInClusters:
                a dictionary of form d['clusterId']{'cell type1', 'cell type2'}

        arguments:
            pandasGroupedData:
                example : df.groupby(['louvain', 'Cell type'])['louvain'].count()

            keys:
                example: grouped.groups.keys()
        '''
        df = anndata.obs.loc[:,['louvain', 'Cell type']]
        grouped = df.groupby(['louvain', 'Cell type'])
        pandasGroupedData=grouped['louvain'].count()
        keys=grouped.groups.keys()

        cellCountsByClusterId = {}
        cellTypesInClusters = {}
        for tup in keys:
            clusterId,cellType = tup
            if cellType not in cellCountsByClusterId :
                cellCountsByClusterId[cellType] = {}

            cellCountsByClusterId[cellType][clusterId] = pandasGroupedData[tup]

            if clusterId not in cellTypesInClusters:
                cellTypesInClusters[clusterId]= set()

            cellTypesInClusters[clusterId].add(cellType)

        return cellCountsByClusterId, cellTypesInClusters


    @staticmethod
    def cellCount(df, cellType):
        rows = (df.loc[:,['Cell type']] == cellType)
        n = rows.sum()
        # n is a pandas series
        return n.values[0]



    @staticmethod
    def annotationProbsForCluster(anndata, cellTypesInClusters, clusterId, cellCountsByClusterId):
        '''
        aedwip: todo:

        arguments:

        returns a list of tuples of following form
            (cellType, p, x, M, n, N)
            p = p-value
            x = number of cells in cluster of cellType
            M = populations size
            n = number of cell types in population
            N = number of cells in the cluster

        '''
        ret = []
        df = anndata.obs.loc[:,['louvain', 'Cell type']]
        M = df.count()[0] # size of population

        cellTypes = cellTypesInClusters[clusterId]
        for cellType in cellTypes:
            # count the number of cells in the cluster
            rows = df['louvain'] == clusterId
            N = sum(rows)

            totalCount = hw2q2.cellCount(df, cellType)
            n = totalCount # number of cell types in populations

            randomVariable = stats.hypergeom(M, n, N)
            x = cellCountsByClusterId[cellType][clusterId]
            pValue = 1.0 - randomVariable.cdf(x)

            ret.append( (cellType, pValue, x, M, n, N) )

        return ret

    @staticmethod
    def bestAnnotation(anndata, cellTypesInClusters, cellCountsByClusterId):
        '''
        uses hypergeometric distribution to best-matching cell type annotation

        returns:
            dataframe with a row for each cluster. the columns are 'cell type' and probablity
        '''
        clusterIds = cellTypesInClusters.keys()
        retDF = pd.DataFrame()

        K_CELLTYPE = 0
        K_P_VALUE = 1

        for clusterId in clusterIds:
            stats = hw2q2.annotationProbsForCluster(anndata, cellTypesInClusters, clusterId, cellCountsByClusterId)

            best = min(stats, key=lambda tup : tup[K_P_VALUE])

            bestDF = pd.DataFrame( data={'clusterId':int(clusterId), 'Cell type':best[K_CELLTYPE],
                                         'p-value':best[K_P_VALUE]}, index=[int(clusterId)] )

            retDF = retDF.append(bestDF)

        return retDF

    @staticmethod
    def getCellsIdxForCluster(anndata, clusterId):
        '''
        returns a list of indices that can be use to select the cells in
        the cluster. the indices are int values that correspond to the
        rows in the numpy array anndata.obsm['X_pca']

        assumes clustering algorithm was run and results stored in
        anndata.obs['louvain']

        arguments:
            anndata
            clusterId: a string
        '''
        numCells = anndata.obs['louvain'].size

        numpyArrayIdx = [i for i in range(numCells)]
        d = {'louvain':anndata.obs['louvain'], 'npIdx':numpyArrayIdx }
        louvainDF = pd.DataFrame(data=d)

        clusterCells = louvainDF['louvain'] == clusterId
        ret = louvainDF.loc[clusterCells, ['npIdx']]
        return ret.values.flatten()

    @staticmethod
    def testGetCellsIdxForCluster(anndata):
        ret = hw2q2.getCellsIdxForCluster(anndata,clusterId='9')
        print("TEST len(ret):{}".format(len(ret)))

    @staticmethod
    def getGeneExpressionSignatureForCluster(anndata, clusterId):
        '''
        returns the centroid for cluster

        assumes:
            1) gene expression values stored in anndata.X

            2) clustering algorithym was run and results stored in
                anndata.obs['louvain']

        arguments:
            anndata
            clusterId: a string
        '''
        cellIndices = hw2q2.getCellsIdxForCluster(anndata ,clusterId)
        # pathways are in gene space not pca(50) space
        # anndata.obsm is a pandas data frame
        # data = anndata.obsm['X_pca'][cellIndices]
        # anndata.X is numpy array
        data = anndata.X[cellIndices, :]

        byColumns = 0
        return np.mean(data, axis=byColumns)

    @staticmethod
    def testGetGeneExpressionSignatureForCluster(anndata):
        ret = hw2q2.getGeneExpressionSignatureForCluster(anndata, clusterId='9')
        # Sample data
        # PCA 50
        #expFirst = np.array([-5.547451, 13.029236, -0.9483415,
        #                    -5.8531246, -2.0250516])
        #expLast = np.array([-0.05719902, -0.16065401, 0.14730875,
        #                    0.00936914, -0.06920523])

        # anndata.XC
        print("TEST ret[0:5]:\n{}".format(ret[0:5]))
        print("TEST ret[-5:]:\n{}".format(ret[-5:]))

    @staticmethod
    def calculateClusterSignatures(anndata):
        '''
        for each cluster calls getGeneExpressionSignatureForCluster()

        returns a dictionary. Key is the cluster id, value is its signature
        '''

        numCells = anndata.obs['louvain'].size
        numpyArrayIdx = [i for i in range(numCells)]
        d = {'louvain':anndata.obs['louvain'], 'npIdx':numpyArrayIdx }
        louvainDF = pd.DataFrame(data=d)

        clusterSigs = {}
        for clusterId in pd.unique(louvainDF['louvain']) :
            sig = hw2q2.getGeneExpressionSignatureForCluster(anndata, clusterId)
            clusterSigs[clusterId] = sig

        return clusterSigs

    @staticmethod
    def rankPathWaysForCluster(anndata, clusterSigs, clusterId):
        '''
        # to speed up testing do not run if results have already been calculated

        returns a tuple (cvsPath, preRes)
            if cvsPath exists preRes == None
        '''
        rnkDf = pd.DataFrame(data={'gene':anndata.var['GeneName-0'],
                                   'score':clusterSigs[clusterId]})

        base = './gseapy.out/prerank_report_GO_Biological_Process_2018_clusterId-'
        path = base + clusterId
        #print("path:{}".format(path))

        # path to output from previous run
        csvPath = path + "/" + 'gseapy.prerank.gene_sets.report.csv'
        #print("cvsPath:{}".format(csvPath))

        # to speed up debugging do not run if results have already been
        # calculated
        ret = None
        exists = os.path.isfile(csvPath)
        if exists:
            # gseapy.gsea.Prerank
            # type(ret):<class 'gseapy.gsea.Prerank'>
            #ret.outdir = csvPath
            ret = (csvPath, None)
        else:
            preRes = gp.prerank(rnk=rnkDf,
                             gene_sets='GO_Biological_Process_2018',
                             processes=4,
                             permutation_num=100, # reduce number to speed up test
                             outdir=path,format='png')
            ret = (csvPath, preRes)

        #print("type(ret):{}".format(type(ret)))
        return ret

    @staticmethod
    def rankPathWays(anndata, clusterSigs, topN=5):
        '''
        returns a panda dataframe with columns Term, nes, cluster id
            Term: pathway
            nes: normalized enrichment score
            cluster id: an integer
        '''
        retDF = pd.DataFrame()
        for clusterId in clusterSigs.keys():
            csvPath, preRes = hw2q2.rankPathWaysForCluster(anndata,
                                             clusterSigs,
                                             clusterId)
            df = pd.read_csv(csvPath)
            df2 = df.sort_values(by='nes', ascending=False)
            df3 = df2.iloc[0:topN, [0, 2]]
            numRows = df3.shape[0]
            cid = [int(clusterId) for j in range(numRows)]
            df3['cluster id'] = cid
            retDF = retDF.append(df3)

        retDF = retDF.sort_values(by='cluster id')
        return retDF




if __name__ == "__main__":
    main()

