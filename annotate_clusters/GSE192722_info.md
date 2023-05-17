We have separate files for each sample, e.g. "GSM5761690_filtered_feature_bc_matrix_Don1.tar.gz". They also provide "filtered_feature_bc_matrix_Pt1" and "filtered_feature_bc_matrix_Pt2". Unzip and extract from tar. They also provide R code on GEO.

``` R
counts <- Read10X_h5(filename = "/Users/nicolelustgarten/Dropbox/Weill Cornell/DJames Lab/Theca Stroma/ATAC Seq/FP1 RNA_ATACSeq Data/filtered_feature_bc_matrix.h5")

data_dir <- "data/counts/GSE192722/filtered_feature_bc_matrix"
list.files(data_dir)
# "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"
data <- Read10X(data.dir = data_dir)
# 10X data contains more than one type and is being returned as a list containing matrices of each type.

seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
seurat_object
# An object of class Seurat
# 36601 features across 5607 samples within 1 assay
# Active assay: RNA (36601 features, 0 variable features)

# "After routine quality control16, 48,147 cells remained for downstream for analysis."
# We only have 5607. Try reading in separately

# There are 8 different samples
data_dirs <- list.files("data/counts/GSE192722/", "(Don|Pt)(?!.*\\.tar$)")
data_dirs <- data_dirs[!grepl(".tar", data_dirs)]

path = "data/counts/GSE192722"

# Loop through files and read into Seurat objects
for(dir in paste(path,data_dirs, sep="/")){
  data <- Read10X(data.dir = dir)
  seurat_object <- CreateSeuratObject(counts = data)
  name <- gsub("data/counts/GSE192722/filtered_feature_bc_matrix_","",dir)
  assign(name, seurat_object)
}



```
