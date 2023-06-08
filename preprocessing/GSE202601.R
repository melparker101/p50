##############################################
# This code is for preparing the snRNA-seq dataset GSE202601
# We need to filter out all of the older donors
# Use conda env 'seurat2h5'
##############################################

# Seurat command list:
# https://satijalab.org/seurat/articles/essential_commands.html

# Load libraries
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(data.table)
library(dplyr)

# Dataset accession
dataset_acc <- "GSE202601"

# Set path
data_path <- paste0("data/counts/",dataset_acc,"/")

# Input files
mat_file <- paste0(data_path,"GSE202601_human_ovary_snRNA-seq_count.rds")
meta_file <- paste0(data_path,"GSE202601_human_ovary_snRNA-seq_metadata.txt")

# Output_files
seurat_out <- paste0(data_path,"counts_GSE202601.h5Seurat")
meta_out <- paste0(data_path,"metadata_GSE202601.csv")

# Load in data
mat <- readRDS(file = mat_file)
meta <- as.data.frame(fread(meta_file, header=T))

### Format metadata ready for input to CELLEX ###
#################################################

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

# Check dimensions
dim(mat)
dim(meta)

# We have 42568 cells for mat and 23534 cells for meta

# Filter mat to only keep the cells from young donors
ids_use <- meta$cell_id
mat <- mat[,ids_use]

# Write metadata as csv
write.csv(meta, meta_out, row.names=F)

### Convert RDS count file to h5ad ###
######################################

# Create a seurat object containing just the raw counts
seurat_obj <- CreateSeuratObject(counts=mat)

# Take a look at the raw counts
GetAssayData(object = seurat_obj, slot = "counts")[1:10,1:10]

# Save counts as a h5ad file
SaveH5Seurat(seurat_obj, filename = seurat_out)
Convert(seurat_out, dest = "h5ad")

# End
