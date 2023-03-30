##############################################################
# Preprocessing GSE213216 scRNA-seq dataset
# melodyjparker14@gmail.com - March 2023
# This code:
# - Reads in the GSE213216 dataset as a Seurat object 
# - Filters to only include cells from unaffected ovaries
# - Converts/saves the filtered Seurat object to a h5ad file
###############################################################

# Load libraries
library(Seurat)

# Dataset accession
dataset_acc <- "GSE213216"

# Set path
data_path <- paste0("counts/",dataset_acc,"/")

# Input files
seurat_in <- paste0(data_path,"_aux.seurat.shared.rds")

# Output files
seurat_out <- paste0(data_path,"counts_","GSE213216",".h5Seurat")

# Load in Seurat object
aux <- readRDS(file = "_aux.seurat.shared.rds")

# Extract metadata
metadata <- aux[[]]

# View colnames of metadata
colnames(metadata)
# metadata[1:3,]

# View cell types
table(metadata$active.cluster)

# View Major class categories
table(metadata$Major.Class)

# Extract raw count data
raw <- GetAssayData(object = aux, slot = "counts")
raw[1:10,1:20]
dim(raw)
# 373,851 cells total, 30354 genes

# Filter for normal ovary
normal_ovary <- subset(aux, subset=Major.Class %in% "Unaffected ovary")

# Extract metadata from filtered the data
ovary_meta <- normal_ovary[[]]

# Check that the metadata only contains data from unaffected ovary
table(ovary_meta$Major.Class)

# Extract raw counts from filtered data
ovary_raw <- GetAssayData(object = normal_ovary, slot = "counts")

# Check the dimension
dim(ovary_raw)
# 30354 51927. 51,927 cells
# ~50,000 cells is what we were expecting from reading the literature

# Save as a h5ad
seurat_obj <- normal_ovary
SaveH5Seurat(seurat_obj, filename = seurat_out)
Convert(seurat_out, dest = "h5ad")
