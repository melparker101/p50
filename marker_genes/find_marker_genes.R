library(Seurat)
library(SeuratDisk)

GSE118127  GSE202601  GSE213216

# Convert("counts/GSE118127/local.h5ad", dest = "h5seurat", overwrite = FALSE)
Convert("counts/GSE118127/local.h5ad", dest = "h5seurat", overwrite = TRUE)
GSE118127 <- LoadH5Seurat("counts/GSE118127/local.h5seurat", verbose = T, misc=F)

LoadH5Seurat("temp.h5seurat", verbose = T)

# https://github.com/mojaveazure/seurat-disk/issues/7

GSE202601

GSE213216

seuratObject <- LoadH5Seurat("example_dir/example_ad.h5Seurat"

# Convert to Seurat ob
Convert("GSE118127/local.h5ad", dest = "h5seurat", overwrite = False)
pbmc3k <- LoadH5Seurat("pbmc3k_final.h5seurat")
pbmc3k
#> An object of class Seurat 
#> 13714 features across 2638 samples within 1 assay 
#> Active assay: RNA (13714 features, 0 variable features)
#>  2 dimensional reductions calculated: pca, umap
