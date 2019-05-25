#!/usr/bin/env python3 -u
# https://stackoverflow.com/a/27534908/4586180
# -u prevents python from buffering stdout. by using -u you can use the 'tee' command

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
import louvain as lv

import logging
from setupLogging import setupLogging
logConfigFile='./test/logging.test.ini.json'
setupLogging( default_path=logConfigFile)

from timeit import default_timer as timer
from datetime import timedelta

def main():
    logger = logging.getLogger(__name__)
    if not os.path.isfile(logConfigFile):
        logger.error("missing logConfigFile file:{}".format(logConfigFile))
        return
        
    anndata = sc.read("PBMC.merged.h5ad")
    
    # run our implementation of nearest neighboors and update anndata
    KnnG(anndata, n_neighbors=12, runPCA=True, nPC=50)
    
    # MacBook Pro (Retina, 15-inch, Late 2013)
    # processor 2.6 GHz Intel Core i7
    # memory 16 GB 1600 MHz DDR
    logger.warning( "BEGIN lv.Louvain.runWithAdata(anndata)")
    start = timer()
    root = lv.Louvain.runWithAdata(anndata)
    end = timer()
    logger.warning("Louvain.runWithAdata execution time:{}"\
                         .format(timedelta(seconds=end-start)))
    logger.warning( "END lv.Louvain.runWithAdata(anndata)\n")
    
    logger.warning("modularity:{}".format(root._Q))
    
    logger.warning("CLUSTER ASSIGMENTS")
    level = root
    while level:
        clusterAssignments = level.getClusterAssigments()
        logger.warning("{}:\n{}".format(level._louvainId, clusterAssignments))
        level = level._leafLouvain

if __name__ == '__main__':
    main()
    
