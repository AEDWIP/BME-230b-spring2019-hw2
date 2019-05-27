Answers to questions:
  hw2q*-nb.pdf                              PDF export of Jupyter notebooks
  BME-230b-Spring-2019-hw2_question*.ipynb  Jupyter notebooks

Python code (each of these have main() methods, so you can run them individually from command-line):
  hw2q*.py                                  Source for functions defined in notebooks
  ARIStatistic.py                           Generates ARI statistics from bbl-knn subsamples
  FStatistic.py                             Calculates F-statistic for cluster evaluation
  bblknn.py                                 Implementation of bbl-knn
  euclid_bbknn.py                           Implementation of bbknn
  euclid_knn.py                             Implementation of knn

Louvain code for 2a:
  louvain.py				    Implementation of phase I and II
  louvainDriverShell.py			    Driver for louvain code

Unit tests:
   testLouvainSimple.py
   testlouvain.py
   testlouvainPhase.py
   testlouvainScanpyIntegration.py
   testDisjointGraphs.py		
   testEuclid_bbknn.py			    * unit test for extra credit	
   testEuclid_knnTest.py		    * unit test for extra credit
   logging.test.ini.json

Required dataset:
  PBMC.merged.h5ad                          PBMC single-cell dataset

Miscellaneous:
  BME-230B-hw2-env.yml                      Conda environment specifications
