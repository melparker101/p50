# Clustering and cell type annotation
For clustering, try to follow the methods from the author's paper. This is usually the Seurat pipeline, which consists of:
-   NormalizeData()
-   FindVariableFeatures()
-   ScaleData()
-   RunPCA(verbose = F)
-   RunUMAP(dims = 1:30)
-   FindNeighbors(reduction = "pca", dims = 1:30)
-   FindClusters(resolution = 0.5)

Use [X. Fan et al.](https://www.nature.com/articles/s41467-019-11036-9) as a reference for annotation.
