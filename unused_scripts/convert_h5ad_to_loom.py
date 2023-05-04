# This code converts h5ad files to loom files
# This was not necessary for this dataset as an RDS file was availible to download
# Python=3.10 works

import anndata
import loompy

# This code converts h5ad to loom
# https://www.biostars.org/p/440922/
adata = anndata.read("counts/GSE118127/local.h5ad")
adata.write_loom("counts/GSE118127/local.h5ad.loom")
