#!/usr/bin/env python
# coding: utf-8

# # BME-230B Spring 2019 HW 2 Question 
# James Casaletto, Andrew Davidson, Yuanqing Xue, Jim Zheng
# 
# - ref
#     * [scanpy.tl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.umap.html)
#     * [scanpy.api.pp.neighbors](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.api.pp.neighbors.html?highlight=neighbors)
#     * [scanpy.pl.umap](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.pl.umap.html#scanpy.pl.umap)
#     * [scanpy.tl.louvain](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.louvain.html#scanpy.tl.louvain)
#     * [GSEAPY: Gene Set Enrichment Analysis in Python. pypi.org](https://pypi.org/project/gseapy/)
#     * [GSEAPY: Gene Set Enrichment Analysis in Python gseapy.readthedocs.io](https://gseapy.readthedocs.io/en/latest/introduction.html)
#     * [anndata](https://anndata.readthedocs.io/en/latest/anndata.AnnData.html)
#         + "uns" stands for unstructured data
#         + "obs" are panda data frame observations 
#         + "obsm key-indexed multi-dimensional observations
#     * [Hypergeometric_distribution](https://en.wikipedia.org/wiki/Hypergeometric_distribution)
#     * [Hypergeometric Tests
# for Gene Lists](http://users.unimi.it/marray/2007/material/day4/Lecture7.pdf)

# In[1]:


from euclid_knn import KnnG
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import os

import pandas as pd
import scanpy.api as sc
import scanpy
print("scanpy.__version__:{}".format(scanpy.__version__))

import scipy.special
import scipy.stats as stats


# ## 2.b. [5 pts] 
# Turn in a UMAP plot of the combined dataset as you did in question #1, but
# this time, color the cells by their Louvain cluster assignments determined for each cell
# within each batch as a different color in each plot.

# In[ ]:


get_ipython().run_cell_magic('time', '', 'anndata = sc.read("PBMC.merged.h5ad")\n\n# run our implementation of nearest neighboors and update anndata\nKnnG(anndata, n_neighbors=12, runPCA=True, nPC=50)')


# In[ ]:


get_ipython().run_cell_magic('time', '', "# running Scanpy's version of Louvian\n# TODO: replace with our implementation\nscanpy.tl.louvain(anndata,\n                  flavor='igraph', \n                  directed=False, \n                  use_weights=True)")


# In[ ]:


plt.figure(figsize=(10,10))
scanpy.tl.umap(anndata)


# In[ ]:


scanpy.pl.umap(anndata, color=['Cell type'])
scanpy.pl.umap(anndata, color=["louvain"])
scanpy.pl.umap(anndata, color=["batch"])


# ## <span style="color:red">AEDWIP TODO</span>
# 
# - after louvain split into  sets 3prime and 5prime
# - create two plots batch vs cell type

# ### 3 and 5 prime cells clusters are different

# In[ ]:


prime3Rows = anndata.obs['Method'] == '10X_3prime'
prime5Rows = anndata.obs['Method'] == '10X_5prime'
print("prime3Rows.sum(): {}".format(prime3Rows.sum()))
print("prime5Rows.shape: {}".format(prime5Rows.sum()))

# save output created by running louvain on combined data set
louvainClustersSave =  anndata.obs.loc[:,['louvain']]
cellTypeSave = anndata.obs.loc[:, ['Cell type']]
batchSave = anndata.obs.loc[:, ['batch']]

# select the 3prime and 5prime data we want to plot
prime3Cluster = louvainClustersSave.loc[prime3Rows,:]
prime5Cluster = louvainClustersSave.loc[prime5Rows,:]

prime3CellType = anndata.obs.loc[prime3Rows,['Cell type']]
prime5CellType = anndata.obs.loc[prime5Rows,['Cell type']]

prime3Batch = anndata.obs.loc[prime3Rows,['batch']]
prime5Batch = anndata.obs.loc[prime5Rows,['batch']]


# In[ ]:


# plot 3 prime
anndata.obs['louvain'] = None # make sure there are not artifacts
anndata.obs.loc[:,['louvain']] = prime3Cluster
scanpy.pl.umap(anndata, color=["louvain"], title="3 prime")

# plot 5 prime
anndata.obs['louvain'] = None # make sure there are not artifacts
anndata.obs.loc[:,['louvain']] = prime5Cluster
scanpy.pl.umap(anndata, color=['louvain'], title="5 prime")


# ### cell typed differ between batchs

# In[ ]:


anndata.obs['louvain'] = None # make sure there are not artifacts
anndata.obs.loc[:,['louvain']] = prime3Cluster
anndata.obs['Cell type'] = None # remove artifacts
anndata.obs.loc[:,['Cell type']] = prime3CellType
scanpy.pl.umap(anndata, color=['Cell type'], title="3 prime")

anndata.obs['louvain'] = None # make sure there are not artifacts
anndata.obs.loc[:,['louvain']] = prime5Cluster
anndata.obs['Cell type'] = None # remove artifacts
anndata.obs.loc[:,['Cell type']] = prime5CellType
scanpy.pl.umap(anndata, color=['Cell type'], title="5 prime")


# ## <span style="color:red">Bug?</span>
# looks like 3 prime has some 5 prime in it. test that 'batch' and 'method' cols agree

# In[ ]:


anndata.obs['louvain'] = None # make sure there are not artifacts
anndata.obs.loc[:,['louvain']] = prime3Cluster
anndata.obs['batch'] = None # remove artifacts
anndata.obs.loc[:,['batch']] = prime3Batch
scanpy.pl.umap(anndata, color=["batch"], title="3 prime")

anndata.obs['louvain'] = None # make sure there are not artifacts
anndata.obs.loc[:,['louvain']] = prime5Cluster
anndata.obs['batch'] = None # remove artifacts
anndata.obs.loc[:,['batch']] = prime5Batch
scanpy.pl.umap(anndata, color=['batch'], title="5 prime")


# In[ ]:


# restore original calculated cluster for combined data set
anndata.obs['louvain'] = None # make sure there are not artifacts
anndata.obs.loc[:,['louvain']] = louvainClustersSave

anndata.obs['Cell type'] = None # make sure there are not artifacts
anndata.obs.loc[:,['Cell type']] = cellTypeSave

anndata.obs['batch'] = None # make sure there are not artifacts
anndata.obs.loc[:,['batch']] = batchSave


# ## 2.c. [5 pts] 
# Turn in a table that lists each cluster and its best-matching cell type
# annotation. The table should contain the cluster number and its best matching cell-type
# annotation based on the hypergeometric analysis.

# In[ ]:


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


# In[ ]:


cellCountsByClusterId, cellTypesInClusters = createCountsDict(anndata)


# In[ ]:


def cellCount(df, cellType):
    rows = (df.loc[:,['Cell type']] == cellType)
    n = rows.sum()
    # n is a pandas series
    return n.values[0]


# In[ ]:


# TODO: AEDWIP: define a class to store results. NO MAGIC NUMBERS
K_CELLTYPE = 0
K_P_VALUE = 1
K_X = 2
K_M = 3
K_n = 4
K_N = 5

def annotationProbsForCluster(anndata, cellTypesInClusters, clusterId):
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
#     print("M:{}".format(M))

    cellTypes = cellTypesInClusters[clusterId]
    for cellType in cellTypes:
        # count the number of cells in the cluster
        rows = df['louvain'] == clusterId
        N = sum(rows)

        totalCount = cellCount(df, cellType)
        n = totalCount # number of cell types in populations

        randomVariable = stats.hypergeom(M, n, N)
        x = cellCountsByClusterId[cellType][clusterId]
        pValue = 1.0 - randomVariable.cdf(x)
        
        ret.append( (cellType, pValue, x, M, n, N) )
        
    return ret
          
ret = annotationProbsForCluster(anndata, cellTypesInClusters, clusterId='0')

# print("\n\n************\nret:\n{}".format(ret))


# In[ ]:


def bestAnnotation(anndata, cellTypesInClusters):
    '''
    uses hypergeometric distribution to best-matching cell type annotation 
    
    returns:
        dataframe with a row for each cluster. the columns are 'cell type' and probablity
    '''
    clusterIds = cellTypesInClusters.keys()
    retDF = pd.DataFrame()
    
    for clusterId in clusterIds:
        stats = annotationProbsForCluster(anndata, cellTypesInClusters, clusterId)
#         print()
#         for s in stats:
#             print("s:\n{}".format(s))

        best = min(stats, key=lambda tup : tup[K_P_VALUE])
#         print("best:\n{}".format(best))

        bestDF = pd.DataFrame( data={'clusterId':int(clusterId),
                                     'Cell type':best[K_CELLTYPE],
                                     'p-value':best[K_P_VALUE]}, 
                                index=[int(clusterId)] )
        retDF = retDF.append(bestDF)
    
    return retDF
        


# In[ ]:


retDF = bestAnnotation(anndata, cellTypesInClusters)  
retDF.sort_values(by=['clusterId'])


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

# In[ ]:


def getCellsIdxForCluster(anndata, clusterId):
    '''
    returns a list of indices that can be use to select the cells in 
    the cluster. the indices are int values that correspond to the
    rows in the numpy array anndata.obsm['X_pca']
    
    assumes clustering algorithym was run and results stored in
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
    
def testGetCellsIdxForCluster(anndata):
    ret = getCellsIdxForCluster(anndata,clusterId='9')
    print("AEDWIP len(ret):{}".format(len(ret)))
    # print("AEDWIP type(ret):{}".format(type(ret)))
    # print ("AEDWIP ret:\n{}".format(ret))
    # count varies slightly from kernal restart to kernal restart
    # sending random.seed() or using , random_state=42 argument to
    # scanpy.louvain() did not resolve this issue
    # assert 747 == len(ret)
    
testGetCellsIdxForCluster(anndata)


# In[ ]:


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
    cellIndices = getCellsIdxForCluster(anndata ,clusterId)
    # pathways are in gene space not pca(50) space
    # anndata.obsm is a pandas data frame
    # data = anndata.obsm['X_pca'][cellIndices]
    # anndata.X is numpy array
    data = anndata.X[cellIndices, :]

    byColumns = 0
    return np.mean(data, axis=byColumns)

def testGetGeneExpressionSignatureForCluster(anndata):
    ret = getGeneExpressionSignatureForCluster(anndata, clusterId='9')
    # PCA 50
    #expFirst = np.array([-5.547451, 13.029236, -0.9483415, 
    #                    -5.8531246, -2.0250516])
    #expLast = np.array([-0.05719902, -0.16065401, 0.14730875, 
    #                    0.00936914, -0.06920523])
    
    # anndata.XC
    print("ret[0:5]:\n{}".format(ret[0:5]))
    print("ret[-5:]:\n{}".format(ret[-5:]))

    # count varies slightly from kernal restart to kernal restart
    # sending random.seed() or using , random_state=42 argument to
    # scanpy.louvain() did not resolve this issue
#     np.testing.assert_array_almost_equal(expFirst, ret[0:5])
#     np.testing.assert_array_almost_equal(expLast, ret[-5:])
    
testGetGeneExpressionSignatureForCluster(anndata)


# In[ ]:


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
        sig = getGeneExpressionSignatureForCluster(anndata, clusterId)
        clusterSigs[clusterId] = sig
    
    return clusterSigs

clusterSigs = calculateClusterSignatures(anndata)


# In[ ]:


def rankPathWaysForCluster(anndata, clusterSigs, clusterId):
    '''
    aedwip
        # to speed up debugging do not run if results have already been 
    # calculated
    
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


# In[ ]:


def rankPathWays(anndata, clusterSigs, topN=5):
    '''
    returns a panda dataframe with columns Term, nes, cluster id
        Term: pathway
        nes: normalized enrichment score
        cluster id: an integer
    '''
    retDF = pd.DataFrame()
    for clusterId in clusterSigs.keys():
        csvPath, preRes = rankPathWaysForCluster(anndata, 
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


# In[ ]:


get_ipython().run_cell_magic('time', '', 'retDF = rankPathWays(anndata, clusterSigs, topN=5)')


# In[ ]:


# https://stackoverflow.com/a/35693013/4586180
# display data frame with out index column
from IPython.display import display, HTML
display(HTML(retDF.to_html(index=False)))

