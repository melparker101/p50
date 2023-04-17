########## Some of this is unnecessary #################

library(Seurat)
library(SeuratDisk)
library(SeuratObject)

GSE118127  GSE202601  GSE213216

# Convert("counts/GSE118127/local.h5ad", dest = "h5seurat", overwrite = FALSE)
Convert("counts/GSE118127/local.gzip.h5ad", dest = "h5seurat", overwrite = FALSE)

GSE118127 <- LoadH5Seurat("counts/GSE118127/local.gzip.h5seurat", verbose = T, misc=F)
GSE118127.loom <- as.loom("counts/GSE118127/local.h5ad.loom", filename = "../output/pbmc3k.loom", verbose = FALSE)

GSE118127 <- Connect(filename = "counts/GSE118127/local.h5ad.loom", mode = "r")
GSE118127   
# GSE118127.seurat <- as.Seurat(GSE118127)
GSE118127.seurat <- as.Seurat(GSE118127, features = "feature_name",cells = "barcode")
GSE118127$close_all()
GSE118127.seurat

# Or we can just directly download the rsd 
# https://cellxgene.cziscience.com/collections/2902f08c-f83c-470e-a541-e463e25e5058

GSE118127.seurat <- readRDS(file = "counts/GSE118127/local.rds")
GSE118127.seurat 

##### We have 19 clusters for GSE118127 #####

# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
GSE118127.seurat <- FindVariableFeatures(GSE118127.seurat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(GSE118127.seurat), 10)
# Plot with and without labels
plot1 <- VariableFeaturePlot(GSE118127.seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# PCA
all.genes <- rownames(GSE118127.seurat)
GSE118127.seurat <- ScaleData(GSE118127.seurat, features = all.genes)
# perform PCA
GSE118127.seurat <- RunPCA(GSE118127.seurat, features = VariableFeatures(object = GSE118127.seurat))
print(GSE118127.seurat[["pca"]], dims = 1:5, nfeatures = 5)

# Visualise PCA
VizDimLoadings(GSE118127.seurat, dims = 1:2, reduction = "pca")
DimPlot(GSE118127.seurat, reduction = "pca")

# UMAP
GSE118127.seurat <- RunUMAP(GSE118127.seurat, dims = 1:10)
DimPlot(GSE118127.seurat, reduction = "umap")



GSE118127.seurat[[c('cluster_id','cell_description','cell_type')]][1:20,]
table(GSE118127.seurat[['cluster_id']])


levels(GSE118127.seurat)

# Set identity classes to an existing column in meta data
Idents(object = GSE118127.seurat) <- "cell_description"
levels(GSE118127.seurat)

out_dir <- "cluster_markers/GSE118127"
dir.create(out_dir)

cell_type_list <- levels(GSE118127.seurat)
for (cell_type in cell_type_list){
  name <- paste(cell_type,"markers",sep="_")
  name <- sub(" ", "_", name)
  name <- sub(")", "", name)
  name <- sub("-", "_", name)
  value <- FindMarkers(GSE118127.seurat, ident.1 = cell_type)
  assign(name, value)
  out <- paste(out_dir,name,sep="/")
  write.table(value,out,sep="\t")
}
  


GSE202601

GSE213216

seuratObject <- LoadH5Seurat("example_dir/example_ad.h5Seurat")
# Convert to Seurat ob
Convert("GSE118127/local.h5ad", dest = "h5seurat", overwrite = False)
pbmc3k <- LoadH5Seurat("pbmc3k_final.h5seurat")
pbmc3k
