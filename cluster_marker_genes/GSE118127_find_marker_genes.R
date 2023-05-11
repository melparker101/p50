##############################################################################
## Find marker genes for dataset GSE118127 (9 clusters and 19 clusters)
## Using a Wilcoxon Rank Sum test
## melodyjparker@gmail.com - Apr 23
##############################################################################

# Directly download the RDS 
# https://cellxgene.cziscience.com/collections/2902f08c-f83c-470e-a541-e463e25e5058

# The object should be filtered and preprocessed already. Their filtering code:
# https://github.com/johnmous/singleCell/blob/master/workflow.Rmd

################################
# Load libraries
################################

library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(dplyr)
library(pals)
library(ggplot2)

################################
# Define functions
################################

# Function for plotting clusters
plotClusters <- function(object,clusters,out){
  n <- length(levels(object))
  clust_no <- paste0(n,"C")
  cluster.cols = as.vector(polychrome(n))
  plot = DimPlot(object = seurat_ob, 
               pt.size = 0.1, 
               raster=FALSE, 
               group.by = clusters,
               cols = cluster.cols,
               label = T)
  ggsave(paste0(out,"/",clusters,"_clusters_",clust_no,".pdf"), width = 10, height = 10)
}

################################
# Set up
################################

# rm(list = ls())

dataset <- "GSE118127"

in_file <- paste0("data/counts/",dataset,"/local.rds")
out_dir <- paste0("cluster_markers/", dataset)

# Create output directory
dir.create(out_dir)

# Load in Seurat object
seurat_ob <- readRDS(file = in_file)
seurat_ob

# Have a look at the metadata
seurat_ob[[c('cluster_id','cell_description','cell_type')]][1:20,]
table(seurat_ob[['cluster_id']])

# cell_description corresponds to cluster_id

# Create output directory
dir.create(out_dir)

# View metadata on feature level
head(seurat_ob[['RNA']][[]])
# Colnames are ensembl ids - we cannot change this
# There is a feature_name column

# The findAllMarkers uses the "data" slot, so change the rownames of this
rownames(seurat_ob@assays$RNA@data) <- seurat_ob[['RNA']][['feature_name']][,1]

# Check that the seurat object has been filtered 
VlnPlot(object = seurat_ob,
        features = c("nGene", "nUMI", "percent.mito"))

# Check that the data is normalised
GetAssayData(object = seurat_ob)[1:10,1:15]

# Check dimensions
dim(GetAssayData(object = seurat_ob))
# 32922 genes

################################
# Find markers
################################

### 1. Use cell_description (9 clusters) ###
############################################

use_col <- "cell_type"

# Set identity classes to the cell description
Idents(object = seurat_ob) <- use_col
levels(seurat_ob)

# Create new directory for results of this cluster set
clust_no <- paste0(length(levels(seurat_ob)),"C")
cluster_dir <- paste(out_dir,clust_no,sep="/")
dir.create(cluster_dir)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
combined_markers <- FindAllMarkers(object = seurat_ob, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)  
View(combined_markers)

# pval = 0 just means the p-value is too small for R to display
# https://github.com/satijalab/seurat/issues/205

# Order the rows by clusters, then by p-values
combined_markers <- combined_markers %>% arrange(as.character(cluster), as.numeric(as.character(p_val)))
combined_markers <- combined_markers %>% relocate(gene) %>% relocate(cluster)
View(combined_markers)

# Write as table
combined_out <- paste0(cluster_dir,"/combined_markers_" , clust_no,".txt")
write.table(combined_markers,combined_out,sep="\t",quote = FALSE)

# Find top 5 markers per cluster by 2-fold change
top5_comb <- combined_markers %>%
        group_by(cluster) %>%
        top_n(n = 5,
              wt = avg_log2FC)
              
# View top 5 markers per cluster
View(top5_comb)

# Write as table
top5_comb_out <- paste0(cluster_dir,"/top5_avglog2FC_comb_markers_",clust_no,".txt")
write.table(top5_comb,top5_comb_out,sep="\t",quote = FALSE)

# Check that there are no clusters with less than 3 cells
table(seurat_ob$cell_type)

# Find markers for each cluster and write to separate tables 
cell_type_list <- levels(seurat_ob)
for (cell_type in cell_type_list){
  # Extract and edit cell type names
  name <- paste(cell_type,"markers",sep="_")
  name <- gsub(" ", "_", name)
  name <- gsub(")", "", name)
  name <- gsub("-", "_", name)
  name <- gsub("/", "_", name)
  # Find markers
  value <- FindMarkers(seurat_ob, ident.1 = cell_type, only.pos = TRUE)
  # Reformat marker results table - order by pval and change col order
  value <- value %>% arrange(as.numeric(as.character(p_val)))
  assign(name, value)
  # Write to file
  out <- paste0(cluster_dir,"/",name,"_",clust_no,".txt")
  write.table(value,out,sep="\t",quote = FALSE)
}

# Plot clusters
# plotClusters(seurat_ob,use_col,cluster_dir)

### 2. Use cell_description (19 clusters) ###
#############################################

use_col <- "cell_description"

# Set identity classes to the cell description
Idents(object = seurat_ob) <- use_col
levels(seurat_ob)

# Create new directory for results of this cluster set
clust_no <- paste0(length(levels(seurat_ob)),"C")
cluster_dir <- paste(out_dir,clust_no,sep="/")
dir.create(cluster_dir)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
combined_markers <- FindAllMarkers(object = seurat_ob, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)  
View(combined_markers)

# Order the rows by clusters, then by p-adjusted values
combined_markers <- combined_markers %>% arrange(as.character(cluster), as.numeric(as.character(p_val)))
combined_markers <- combined_markers %>% relocate(gene) %>% relocate(cluster)
View(combined_markers)

# Write as table
combined_out <- paste0(cluster_dir,"/combined_markers_" , clust_no,".txt")
write.table(combined_markers,combined_out,sep="\t",quote = FALSE)

# Find top 5 markers per cluster by 2-fold change
top5_comb <- combined_markers %>%
        group_by(cluster) %>%
        top_n(n = 5,
              wt = avg_log2FC)
View(top5_comb)

# Write as table
top5_comb_out <- paste0(cluster_dir,"/top5_avglog2FC_comb_markers_",clust_no,".txt")
write.table(top5_comb,top5_comb_out,sep="\t",quote = FALSE)

# Find markers for each cluster and write to separate tables 
cell_type_list <- levels(seurat_ob)
for (cell_type in cell_type_list){
  # Extract and edit cell type names
  name <- paste(cell_type,"markers",sep="_")
  name <- gsub(" ", "_", name)
  name <- gsub(")", "", name)
  name <- gsub("-", "_", name)
  name <- gsub("/", "_", name)
  value <- FindMarkers(seurat_ob, ident.1 = cell_type, only.pos = TRUE)
  # Reformat marker results table - order by pval and change col order
  value <- value %>% arrange(as.numeric(as.character(p_val)))
  assign(name, value)
  # Write to file
  out <- paste0(cluster_dir,"/",name,"_",clust_no,".txt")
  write.table(value,out,sep="\t",quote = FALSE)
}

# Plot clusters
# plotClusters(seurat_ob,use_col,cluster_dir)

################################
# End
################################
