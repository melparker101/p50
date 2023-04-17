## This code removes 
import anndata

adata = anndata.read("counts/GSE118127/local.h5ad")
adata.write("counts/GSE118127/local.gzip.h5ad", compression="gzip")
