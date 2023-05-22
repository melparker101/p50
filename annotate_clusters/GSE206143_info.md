We have a matrix, barcode and feature file for each of the 4 samples. Download and extract from tar.

``` bash
mkdir data_dir_IVF1 data_dir_IVF2 data_dir_IVF3 data_dir_IVF4
for f in *IVF1*.gz ; do mv $f "data_dir_IVF1/"$(echo "$f" | sed 's/.*_//'); done
for f in *IVF2*.gz ; do mv $f "data_dir_IVF2/"$(echo "$f" | sed 's/.*_//'); done
for f in *IVF3*.gz ; do mv $f "data_dir_IVF3/"$(echo "$f" | sed 's/.*_//'); done
for f in *IVF4*.gz ; do mv $f "data_dir_IVF4/"$(echo "$f" | sed 's/.*_//'); done
```
 Organise like this:
|-- data_dir_IVF1
|   |-- barcodes.tsv
|   |-- features.tsv
|   `-- matrix.mtx
|-- data_dir_IVF2
|   |-- barcodes.tsv
|   |-- features.tsv
|   `-- matrix.mtx
|-- data_dir_IVF3
|   |-- barcodes.tsv
|   |-- features.tsv
|   `-- matrix.mtx
`-- data_dir_IVF4
    |-- barcodes.tsv
    |-- features.tsv
    `-- matrix.mtx

"A UMAP shows 19 clusters (0–18) of 7609 cells via unsupervised clustering analyses."

``` R
# Load libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)

# Define dataset, data path and seurat outfile
dataset <- "GSE206143"
path <- paste0("data/counts/",dataset)
out_seurat <- paste0(path, "/", dataset, "_counts.h5Seurat")

# Make a list of the data directories 
data_dirs <- list.files(path, "data_dir")

# Make an empty vector to add sample names to
samples <- c()

# Loop through samples
for(dir in data_dirs){

  # Pick sample name
  name <- gsub("data_dir_", "", dir)
  
  # Read in data
  data <- Read10X(data.dir = paste(path,dir,sep="/"))  
  
  # Create Seurat object
  seurat_object <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
  
  # Process Seurat object
  seurat_object <- subset(x = seurat_object, subset = nFeature_RNA > 500)
  seurat_object <- NormalizeData(object = seurat_object, verbose = FALSE)
  seurat_object <- FindVariableFeatures(object = seurat_object, selection.method = "vst", nfeatures = 2000)
  
  # Assign sample name to seurat object
  assign(name, seurat_object)
  samples <- append(samples,name)
  
  # Clean up
  rm(data)
  rm(seurat_object)
}

# Create a command for merging seurat objects
# Create a list of P2-P5 separated by commas
merge_list <- paste(samples[2:length(samples)], collapse = ",")

# "add.cell.ids = patient_id" to avoid duplicate cell names
mergeSamples <- paste0("merged_ob <- merge(",samples[1],", y = c(", merge_list, "), add.cell.ids = samples, project = \"ovaries, min.cells=3, min.features = 200 \")")

# Run command to merge Seurat objects
merged_ob <- eval(parse(text = mergeSamples))
merged_ob
# 8318 cells

# Extract raw data
raw_counts <- GetAssayData(object = merged_ob, slot = "counts")

## Find Anchors
anchors <- FindIntegrationAnchors(object.list = list(IVF1, IVF2, IVF3, IVF4), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = combined) <- "integrated"

## visualization and clustering
combined <- ScaleData(object = combined, verbose = FALSE)
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
DefaultAssay(object = combined) <- "RNA"

# saveRDS(combined, file="integrated_IVFs.rds")


# readRDS("integrated_IVFs.rds")

p1 <- DimPlot(object = combined, reduction = "umap", group.by = "patient")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1)
plot_grid(p2)
DimPlot(combined, reduction = "umap", pt.size = 1, label = TRUE, label.size = 9) + NoLegend()

## UMAPs showing 4 individual IVF patients
DimPlot(object = combined, reduction = "umap", split.by = "patient")


#Cells showing a minimum expression of 500 genes were filtered, and 
#genes showing the expression in a minimum of five cells were applied.  # it says 3 in the code
#After log-normalization using NormalizeData(), 
#four datasets were integrated following FindIntegrationAnchors() and IntegrateData() 
#procedures with the first 20 principal components for weighting. Filtered cells were dimensionally reduced and visualized on a Uniform Manifold Approximation and Projection (UMAP).25 Unsupervised clustering analyses on the UMAP embedding were conducted using #FindClusters() with “resolution = 0.05,” and cells were divided by four individual IVF patient samples in the same UMAP.




#################################






```
