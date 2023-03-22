# This code is for preparing the snRNA-seq dataset GSE202601

# Load libraries
library(Seurat)
library(data.table)
library(dplyr)

# Input files
mat_file <- "GSE202601/GSE202601_human_ovary_snRNA-seq_count.rds"
meta_file <- "GSE202601/GSE202601_human_ovary_snRNA-seq_metadata.txt"
meta_out <- "metadata_GSE202601.csv"

# Load in data
mat <- readRDS(file = mat_file)
meta <- as.data.frame(fread(meta_file, header=T))

# Format metadata ready for input to CELLEX
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
#dim(mat)
dim(meta)

# We have 42568 cells for both

# Write metadata as csv
write.csv(meta, meta_out)
