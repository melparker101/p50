##############################################################
# Preprocessing GSE213216 scRNA-seq dataset
# melodyjparker14@gmail.com - March 2023
# Use at least three cores, preferably more
# $PWD = p50
# This code:
# - Reads in the GSE213216 dataset as a Seurat object 
# - Filters to only include cells from unaffected ovaries
# - Converts/saves the filtered Seurat object to a h5ad file
##############################################################

# Load libraries
library(Seurat)
library(SeuratDisk)
library(SeuratData)

# Dataset accession
dataset_acc <- "GSE213216"

# Set path
data_path <- paste0("data/counts/",dataset_acc,"/")

# Input files
seurat_in <- paste0(data_path,"_aux.seurat.shared.rds")

# Output files
seurat_out <- paste0(data_path,"counts_",dataset_acc,".h5Seurat")

# Load in Seurat object
aux <- readRDS(file = seurat_in)

# Have a look at the colnames of the metadata
colnames(aux[[]])

# Create a new metadata column and delete the old column to replace . with _ for our columns of interest
aux$major_class <- aux$Major.Class
aux$Major.Class <- NULL

aux$active_cluster <- aux$active.cluster
aux$active.cluster <- NULL

colnames(aux[[]])

# View how many cells we have of each cell type
table(aux$active_cluster)

# View how many cells we have from each sample type category
table(aux$major_class)

# Filter for normal ovary samples
ovary <- subset(aux, subset=major_class %in% "Unaffected ovary")

dim(ovary)
# 30354 51927
# 51,927 cells is what we were expecting for unaffected ovary from the literature (Fig. 2)

# Check that the metadata only contains data from unaffected ovary
unique(ovary$major_class)

# Now filter for pre-menopausal patients
# Take a look at the patients
table(ovary$Patient.No.)
# Patient IDs are 6, 14, 15, 18

# From Supplementary Table two we can see that patients 14 and 15 are pre-menopausal
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01254-1/MediaObjects/41588_2022_1254_MOESM4_ESM.xlsx
ovary <- subset(ovary, subset=Patient.No. %in% c(14,15))

# Check that it has filtered properly
unique(ovary$Patient.No.)
# 14 15

dim(ovary)
# 30354 22219
# We now have 22,219 cells

# Change metadata to characters so that it doesn't convert to h5ad as integers
i <- sapply(ovary@meta.data, is.factor)
ovary@meta.data[i] <- lapply(ovary@meta.data[i], as.character)

##############################################################

# Extract raw counts from filtered data
raw <- GetAssayData(object = ovary, slot = "counts")

# Extract metadata from filtered the data
metadata <- ovary[[]]

# Create a new Seurat object using only the raw data - otherwise raw count data won't convert to anndata properly
seurat_obj <- CreateSeuratObject(counts=raw, metadata=meta.data)

# Save as a h5ad
SaveH5Seurat(seurat_obj, filename = seurat_out)
Convert(seurat_out, dest = "h5ad")

# Done
