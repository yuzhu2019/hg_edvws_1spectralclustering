# hg_edvws_1spectralclustering

The code for the paper "Hypergraphs with Edge-Dependent Vertex Weights: p-Laplacians and Spectral Clustering" by Yu Zhu and Santiago Segarra.

The proposed algorithm (Algorithm 2 in our paper) is implemented in the following two files (under the folder 'matlab'):

(1) reducible_hypergraph_partition.cpp       --  The inner problem is solved via FISTA (see Algorithm 3 in our paper)

(2) reducible_hypergraph_partition_pdhg.cpp  --  The inner problem is solved via PDHG  (see Algorithm 4 in our paper)
