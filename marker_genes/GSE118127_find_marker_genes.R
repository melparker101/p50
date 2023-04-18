# We have 19 clusters for GSE118127
# Find markers for these and write output as tables
library(Seurat)
library(SeuratDisk)
library(SeuratObject)

in_file <- "counts/GSE118127/local.rds"
out_dir <- "cluster_markers/GSE118127"

# Directly download the RDS 
# https://cellxgene.cziscience.com/collections/2902f08c-f83c-470e-a541-e463e25e5058

# Load in Seurat object
GSE118127.seurat <- readRDS(file = in_file)
GSE118127.seurat 

# Have a look at the metadata
GSE118127.seurat[[c('cluster_id','cell_description','cell_type')]][1:20,]
table(GSE118127.seurat[['cluster_id']])

# Check current identity classes
levels(GSE118127.seurat)

# Set identity classes to the cell description column in meta data
Idents(object = GSE118127.seurat) <- "cell_description"
levels(GSE118127.seurat)

# Create output directory
dir.create(out_dir)

# Find markers and write tables for each cluster
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
