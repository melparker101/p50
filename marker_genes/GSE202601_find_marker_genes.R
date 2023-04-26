##############################################################################
## Find marker genes for dataset GSE202601 (8 clusters)
## Using a Wilcoxon Rank Sum test
## melodyjparker@gmail.com - Apr 23
##############################################################################

# Use the original rds file and metadata
# The RDS contains the raw counts, so this needs preprocessing/basic filtering

# https://github.com/satijalab/seurat/issues/678
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# Their data processing code:
# https://github.com/ChenJin2020/The-regulatory-landscapes-of-human-ovarian-ageing/blob/main/Data_Processing.R

################################
# Load libraries
################################
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(pals)
library(data.table)

################################
# Set up 
################################

# rm(list = ls())

# Dataset
dataset <- "GSE202601"

# Data path
data_path <- paste0("counts/",dataset,"/")

out_dir <- paste0("cluster_markers/", dataset)

# Create output directory
dir.create(out_dir)

# Input files
mat_file <- paste0(data_path,"GSE202601_human_ovary_snRNA-seq_count.rds")
meta_file <- paste0(data_path,"GSE202601_human_ovary_snRNA-seq_metadata.txt")

# Load in data
mat <- readRDS(file = mat_file)
meta <- as.data.frame(fread(meta_file, header=T))

################################
# Filtering
################################

# Filter out old donors
meta = meta[meta$group == "young",]

# Check the ages of the remaining donors
table(meta$age)

# See 'Metadata format' in https://github.com/perslab/CELLEX
meta <- meta[,c(1,5)]
colnames(meta) <- c("cell_id","abbr")

# Make a table to map cell type abbreviations to full names
abbr <- names(table(meta$abbr))
# "BEC"  "EpiC" "GC"   "IC"   "LEC"  "SC"   "SMC"  "TC"
cell_type <- c("blood_vessel_endothelial_cells", "epithelial_cells", "granulosa_cells", "immune_cells", "lymphatic_endothelial_cells", "stromal_cells", "smooth_muscle_cells","theca_cells")
map <- data.frame(abbr, cell_type)

# Replace cell type abbreviations with their full names
meta <- inner_join(meta, map, by=c("abbr"))[,c(1,3)]

# Filter mat to only keep the cells from young donors
ids_use <- meta$cell_id
mat <- mat[,ids_use]

# Check dimensions
dim(mat)
dim(meta)

################################
# Set up
################################

# Use the original code from the dataset to process and cluster and also use the Seurat tutorial
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# https://github.com/ChenJin2020/The-regulatory-landscapes-of-human-ovarian-ageing/blob/main/Data_Processing.R

# Format metadata
rownames(meta) <- meta$cell_id
meta <- meta[,2,drop=FALSE]

# Create Seurat object
seurat_ob <- CreateSeuratObject(counts = mat, project = "ovary", min.cells = 3, min.features = 200, meta.data = meta)

# calculate mitochondia percentage of cell
seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")

# Plot mitochondria percent and nfeatures
plot1 <- FeatureScatter(seurat_ob, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter data (though it looks like it is already filtered)
seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

# The counts matrix we created the seurat object with was raw counts
mat[1:20,1:20]

# Normalise
seurat_ob <- NormalizeData(seurat_ob, normalization.method = "LogNormalize", scale.factor = 10000)

# View normalised counts
GetAssayData(seurat_ob)[1:20,1:20]

################################
# Find markers
################################

### 1. Use active clusters (8 clusters) ###
###########################################

# Set indentity classes as the active clusters
use_col <- "cell_type"
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

# Order the rows by clusters, then by p-values
combined_markers <- combined_markers %>% arrange(as.character(cluster), as.numeric(as.character(p_val)))
combined_markers <- combined_markers %>% relocate(gene) %>% relocate(cluster)
View(combined_markers)

# Write as table
combined_out <- paste0(cluster_dir,"/combined_markers_" , clust_no,".txt")
write.table(combined_markers,combined_out,sep="\t",quote = FALSE)

# Find top 5 markers per cluster by log2 fold change
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
  # Extract cell type names
  name <- paste(cell_type,"markers",sep="_")
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
