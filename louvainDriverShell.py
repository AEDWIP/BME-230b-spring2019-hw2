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

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    ref: https://docs.python.org/3/library/argparse.html

    '''
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        
        Implements a parser to interpret the command line argv string using argparse.
        
        ref: restriced float https://stackoverflow.com/a/12117065/4586180
        '''
        
        import argparse
        
        def str2bool(v):
            """
            https://stackoverflow.com/a/43357954/4586180
            """
            if v.lower() in ('yes', 'true', 't', 'y', '1'):
                return True
            elif v.lower() in ('no', 'false', 'f', 'n', '0'):
                return False
            else:
                raise argparse.ArgumentTypeError('Boolean value expected.')
                
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output' 
                                             )
      
        self.parser.add_argument('-c', '--clusterAssigmentOutputFile', default=None)
        
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)
  

class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg 
        
        
########################################################################
# CommandLine
########################################################################
def main(myCommandLine=None):
    logger = logging.getLogger(__name__)
    if not os.path.isfile(logConfigFile):
        logger.error("missing logConfigFile file:{}".format(logConfigFile))
        return
        
    logger.warning("BEGIN")
    
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else :
        myCommandLine = CommandLine(myCommandLine) # interpret the list passed from the caller of main as the commandline.

    try:
        print("command line args:{}".format(myCommandLine.args))

#         if myCommandLine.args.requiredBool:
#             print ('requiredBool is', str(myCommandLine.args.requiredBool) )  ## this is just an example
#         else :
#             pass
#         raise Usage ('testing')  # this is an example of how to raise a Usage exception and you can include some text that will get printed. Delete this is you dont need it

    except Usage as err:
        print (err.msg)
         

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
    

    logger.warning("clustering completed successfully root level modularity:{}".format(root._Q))
    
    clusterAssigmentOutputFile = myCommandLine.args.clusterAssigmentOutputFile 
    if not clusterAssigmentOutputFile:
        clusterAssigmentOutputFile = "./louvainDriverShell.py.out"
        logger.warn("clusterAssigment output will be written to:{}".format(clusterAssigmentOutputFile))
        
    logger.info("clusterAssigmentOutputFile:{}".format(clusterAssigmentOutputFile))
            
    with open(clusterAssigmentOutputFile, "w") as f:
        level = root
        while level:
            clusterAssignments = level.getClusterAssigments()
            numCluster = level.countClusters()
            msg = "clustering completed successfully cluster assignments for"
            hdrMsg = "######## {} level:{} numClusters:{}:"\
                .format(msg, level._louvainId, numCluster)
            f.write(hdrMsg)
            f.write(clusterAssignments)
            level = level._leafLouvain
        
    logger.warning("END")
        

if __name__ == '__main__':
    main()
    
