##############################################################################
## Find marker genes for dataset GSE213216 (9 clusters)
## Using a Wilcoxon Rank Sum test
## melodyjparker@gmail.com - Apr 23
##############################################################################

# Use the original _aux.seurat.shared.rds because it contains the normalised counts
# They used SCTransform to normalise the data
# Check if this is still normally distributed for our filtered subset
# Only use active clusters (there are too many seurat clusters)

# Run this for
# 1. The already processed (SCT-normalized) dataset. Filter and then find markers
# 2. The filtered (for CELLEX) raw counts. Normalise these (log normalization), then find markers

################################
# Load libraries
################################
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(pals)

################################
# Functions
################################
# Function for plotting clusters
plotClusters <- function(object,clusters,out){
  n <- length(levels(object))
  clust_no <- paste0(n,"C")
  # cluster.cols = as.vector(as.vector(tol.rainbow(n)))
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
# 1. SCTransform normalisation
################################
##### Set up #####
# rm(list = ls())

# Dataset accession
dataset <- "GSE213216"

# Out dir
out_dir <- paste0("cluster_markers/", dataset)

# Create output directory
dir.create(out_dir)

# Input file - dataset is processed (already normalised) and unfiltered (for the samples we want to keep)
seurat_in <- paste0("counts/",dataset,"/_aux.seurat.shared.rds")

# Load in Seurat object
seurat_ob <- readRDS(file = seurat_in)
seurat_ob

# View the total RNA count per cell and group by clusters
# VlnPlot(object = seurat_ob, features = c("nFeature_RNA"), group.by = c('seurat_clusters'))
VlnPlot(object = seurat_ob, features = c("nFeature_RNA"), group.by = c('active.cluster'))

##### Filtering #####

# Have a look at the colnames of the metadata
colnames(seurat_ob[[]])

# Create a new metadata column and delete the old column to replace . with _ for our columns of interest
seurat_ob$major_class <- seurat_ob$Major.Class
seurat_ob$Major.Class <- NULL

seurat_ob$active_cluster <- seurat_ob$active.cluster
seurat_ob$active.cluster <- NULL

colnames(seurat_ob[[]])

# Filter for normal ovary samples
seurat_ob <- subset(seurat_ob, subset=major_class %in% "Unaffected ovary")

dim(seurat_ob)
# 30354 51927
# 51,927 cells is what we were expecting for unaffected ovary from the literature (Fig. 2)

# Check that the metadata only contains data from unaffected ovary
unique(seurat_ob$major_class)

# Now filter for pre-menopausal patients
# Take a look at the patients
table(seurat_ob$Patient.No.)
# Patient IDs are 6, 14, 15, 18

# From Supplementary Table two we can see that patients 14 and 15 are pre-menopausal
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01254-1/MediaObjects/41588_2022_1254_MOESM4_ESM.xlsx
seurat_ob <- subset(seurat_ob, subset=Patient.No. %in% c(14,15))

# Check that it has filtered properly
unique(seurat_ob$Patient.No.)
# 14 15

dim(seurat_ob)
# 30354 22219
# We now have 22,219 cells

# See how many clusters there are for different cluster types
dim(unique(seurat_ob[['seurat_clusters']]))[1]
seurat_ob[[c('seurat_clusters','active_cluster')]][1:20,]
dim(unique(seurat_ob[['active_cluster']]))[1]

##### Find markers #####

# Use active clusters (9 clusters)
# Cluster names are character
# Number of clusters = 9
norm_method = "SCT"

# Set indentity classes as the active clusters
use_col <- "active_cluster"
Idents(object = seurat_ob) <- use_col
levels(seurat_ob)

# Create new directory for results of this cluster set
clust_no <- paste0(length(levels(seurat_ob)),"C_",norm_method)
cluster_dir <- paste(out_dir,clust_no,sep="/")
dir.create(cluster_dir)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
combined_markers <- FindAllMarkers(object = seurat_ob, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)  
# View(combined_markers)

# Order the rows by clusters, then by p-val
combined_markers <- combined_markers %>% arrange(as.character(cluster), as.numeric(as.character(p_val)))
combined_markers <- combined_markers %>% relocate(gene) %>% relocate(cluster)
# View(combined_markers)

# Write as table
combined_out <- paste0(cluster_dir,"/combined_markers_" , clust_no,".txt")
write.table(combined_markers,combined_out,sep="\t",quote = FALSE)

# Find top 5 markers per cluster
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

################################
# Log normalisation
################################
##### Set up #####

rm(list = ls())

# Dataset accession
dataset <- "GSE213216"

# Out dir
out_dir <- paste0("cluster_markers/", dataset)

# Input file - raw, filtered dataset (created for input to CELLEX)
seurat_in <- paste0("counts/", dataset, "/counts_GSE213216.h5Seurat")

# Read in Seurat object (this is already filtered for the samples we want to use)
seurat_ob <- LoadH5Seurat(seurat_in)
seurat_ob

# Check that this is already filtered
dim(seurat_ob)
# 30354 22219

# View raw counts
GetAssayData(seurat_ob)[1:20,1:20]

# Normalise
seurat_ob <- NormalizeData(seurat_ob, normalization.method = "LogNormalize", scale.factor = 10000)

# View normalised counts
GetAssayData(seurat_ob)[1:20,1:20]

##### Find markers #####

# Use active clusters (9 clusters)
# Cluster names are character
# Number of clusters = 9
norm_method = "log"

# Set indentity classes as the active clusters
use_col <- "active_cluster"
Idents(object = seurat_ob) <- use_col
levels(seurat_ob)

# Create new directory for results of this cluster set
clust_no <- paste0(length(levels(seurat_ob)),"C_",norm_method)
cluster_dir <- paste(out_dir,clust_no,sep="/")
dir.create(cluster_dir)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
combined_markers <- FindAllMarkers(object = seurat_ob, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)  
# View(combined_markers)

# Order the rows by clusters, then by p-val
combined_markers <- combined_markers %>% arrange(as.character(cluster), as.numeric(as.character(p_val)))
combined_markers <- combined_markers %>% relocate(gene) %>% relocate(cluster)
# View(combined_markers)

# Write as table
combined_out <- paste0(cluster_dir,"/combined_markers_" , clust_no,".txt")
write.table(combined_markers,combined_out,sep="\t",quote = FALSE)

# Find top 5 markers per cluster
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
  # Find markers
  value <- FindMarkers(seurat_ob, ident.1 = cell_type, only.pos = TRUE)
  # Reformat marker results table - order by pval and change col order
  value <- value %>% arrange(as.numeric(as.character(p_val)))
  assign(name, value)
  # Write to file
  out <- paste0(cluster_dir,"/",name,"_",clust_no,".txt")
  write.table(value,out,sep="\t",quote = FALSE)
}

################################
# End
################################
