## This code removes LZF filter
# This was not the issue for this dataset

import anndata

# https://github.com/mojaveazure/seurat-disk/issues/7
adata = anndata.read("counts/GSE118127/local.h5ad")
adata.write("counts/GSE118127/local.gzip.h5ad", compression="gzip")

