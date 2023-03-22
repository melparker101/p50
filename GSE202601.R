# This code is for preparing the snRNA-seq dataset GSE202601
####### TIDY THIS UP ########

# Seurat command list:
# https://satijalab.org/seurat/articles/essential_commands.html

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




# Use this CELLEX tutorial to generate 
file from rsd
# https://github.com/perslab/CELLEX/blob/master/tutorials/demo_moca_100k.ipynb

library(Seurat)
library(tidyverse)
library(here)
library(data.table)
library(RLinuxModules)
library(loomR)

cds <- mat

### Get cds cell annotations 
### These annotations does not contain 
### all the information in cell_annotate.csv.gz (e.g. Sub_trajectory_Pseudotime is missing)
df.cell = Biobase::pData(cds) # cell meta-data | same as cds@phenoData@data
df.cell <- df.cell %>% select(sample) %>% as_tibble()

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts=mat, meta.data=meta)

m <- as(mat, "dgTMatrix")
# as(<dgCMatrix>, "dgTMatrix") is deprecated since Matrix 1.5-0; do as(., "TsparseMatrix") instead
m2 <- as(mat, "TsparseMatrix")

# View counts
seurat_obj@assays$RNA@counts[1:10,1:20]

library(SeuratDisk)
library(loomR)

loom <- as.loom(seurat_obj, filename = "GSE202601.loom", verbose = FALSE)


writeMMgz <- function(x, file) {
  mtype <- "real"
  if (is(x, "ngCMatrix")) {
    mtype <- "integer"
  }
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mtype),
      sprintf("%s %s %s", x@Dim[1], x@Dim[2], length(x@x))
    ),
    gzfile(file)
  )
  data.table::fwrite(
    x = summary(x),
    file = file,
    append = TRUE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}

library(Matrix)

writeMMgz(m, "GSE202601/GSE202601_counts.mtx.gz")
m_header <- colnames(m)
write.table(m_header,"GSE202601/matrix_col_names.txt",row.names=FALSE,col.names=FALSE,sep="\t", quote = FALSE)


### First try rsd to seurat to Hdf5
rownames(meta) <- meta[,1]
meta <- meta[-1]
seurat_obj <- CreateSeuratObject(counts=mat, meta.data=meta)
head(seurat_obj@meta.data)

library(SeuratData)

SaveH5Seurat(pbmc3k.final, filename = "GSE202601.h5Seurat")
Convert("GSE202601.h5Seurat", dest = "h5ad")


# python code
import h5py
filename = matrix_GSE202601.h5

with h5py.File(filename, "r") as f:
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
    print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    a_group_key = list(f.keys())[0]

h5 = h5py.File(filename,'r')

futures_data = h5['futures_data']  # VSTOXX futures data
options_data = h5['options_data']  # VSTOXX call option data

h5.close()

import numpy as np
import pandas as pd
df = pd.DataFrame(np.array(h5py.File(path)['variable_1']))

>>> hf = h5py.File('matrix_GSE202601.h5', 'r')
>>> data = hf.get('dataset_name').value # `data` is now an ndarray.
data = hf.get('dataset_name').value # `data` is now an ndarray.

