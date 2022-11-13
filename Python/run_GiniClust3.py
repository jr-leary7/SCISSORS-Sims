def run_gc3(count_mat, res_val):
  # import libraries 
  import anndata           # annotated single cell data
  import numpy as np       # matrix algebra
  import pandas as pd      # DataFrames
  import scanpy as sc      # ScanPy
  import giniclust3 as gc  # GiniClust3
  # create AnnData object
  adata = anndata.AnnData(X=count_mat.T)
  sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
  # calculate clustering 
  gc.gini.calGini(adata, selection='p_value', p_value=0.001)
  gc.fano.calFano(adata, method='scanpy')
  adataGini = gc.gini.clusterGini(adata, resolution=res_val, neighbors=5)
  adataFano = gc.fano.clusterFano(adata, resolution=res_val, neighbors=15)
  consensusCluster = {}
  consensusCluster['giniCluster'] = np.array(adata.obs['rare'].values.tolist())
  consensusCluster['fanoCluster'] = np.array(adata.obs['fano'].values.tolist())
  gc.consensus.generateMtilde(consensusCluster)
  gc.consensus.clusterMtilde(consensusCluster)
  return consensusCluster 
