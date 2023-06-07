We have a matrix, barcode and feature file for each of the 4 samples. Download and extract from tar.

``` bash
mkdir data_dir_IVF1 data_dir_IVF2 data_dir_IVF3 data_dir_IVF4
for f in *IVF1*.gz ; do mv $f "data_dir_IVF1/"$(echo "$f" | sed 's/.*_//'); done
for f in *IVF2*.gz ; do mv $f "data_dir_IVF2/"$(echo "$f" | sed 's/.*_//'); done
for f in *IVF3*.gz ; do mv $f "data_dir_IVF3/"$(echo "$f" | sed 's/.*_//'); done
for f in *IVF4*.gz ; do mv $f "data_dir_IVF4/"$(echo "$f" | sed 's/.*_//'); done

Organise like this:
```text
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
```
"A UMAP shows 19 clusters (0–18) of 7609 cells via unsupervised clustering analyses."

``` R
# Code copied and edited from https://github.com/nurungji82/scRNA-seq_of_IVF_samples/blob/main/R%20script%20for%204%20human%20follicular%20aspirates.R

# Load libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)

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

## Find Anchors
anchors <- FindIntegrationAnchors(object.list = list(IVF1, IVF2, IVF3, IVF4), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = combined) <- "integrated"

## Visualisation and clustering
combined <- ScaleData(object = combined, verbose = FALSE)
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
DefaultAssay(object = combined) <- "RNA"

DimPlot(combined, reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 5) + NoLegend()

# saveRDS(combined, file="integrated_IVFs.rds")
# readRDS("integrated_IVFs.rds")

DimPlot(combined, reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 5) + NoLegend()

## UMAPs showing 4 individual IVF patients
# DimPlot(object = combined, reduction = "umap", split.by = "patient")

# Plot by sample/patient
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

length(table(Idents(combined)))

# After filtering we have the same number of cells (7609)
# After clustering we do not have the same number of cells... we have 18. It's okay, just carry on

###########################################################

## Cell-specific markers for re-clustering
###Follicular cells vs. leukocytes vs. RBCs
FeaturePlot(combined, features = c("CYP11A1", "PTPRC", "HBA1"))

############################

###Follicular cells
VlnPlot(combined, features = c("CYP11A1", "HSD3B2"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CYP11A1", "HSD3B2"))
# Expressed in clusters 0,2,10,11,12 and 17

###Granulosa cells
VlnPlot(combined, features = c("CYP19A1", "INHA"), slot = "counts", log = TRUE)
# From these plots we can predict that clusters 0,2,10,11,12 are GC (0,2,9,11,12 in theirs)
# CCL20 and CXCL12
VlnPlot(combined, features = c("CCL20", "CXCL12"), slot = "counts", log = TRUE)

###Theca cells
VlnPlot(combined, features = c("COL1A1", "COL1A2", "COL3A1"), slot = "counts", log = TRUE)
# 17 is theca (18 in theirs)

############################

###Leukocytes, RBCs, & endothelial cells
VlnPlot(combined, features = c("PTPRC", "HBA1", "TIE1", "VWF"), slot = "counts", log = TRUE)
# Leukocytes PTPRC: 2,3,4,5,6,7,8,9,13,14,15,16
# RBC HBA1: 15
# Endothelial cells: TIE1, VWF: none

###Subclustering leukocytes
####Cytotoxic & Helper T cells
VlnPlot(combined, features = c("CD3E", "CD8A", "CD4"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CD3E", "CD8A", "CD4"))
# Just go with cluster 1 for helper T - it matches their diagram and CD3E and CD4 (cluster 1 for them)
# Again, go for 7 for cytotoxin because it matches their diagram and (4 for them) CD3E and CD8A

####NK cells
VlnPlot(combined, features = c("NKG7", "NCR1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("NKG7", "NCR1"))
# Cluster 6 (6 for theirs)

####NKT cells
VlnPlot(combined, features = c("NKG7", "CD3E","CD8A"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("NKG7", "CD3E","CD8A"))
# Cluster 5 (7 for theirs)

####B cells
VlnPlot(combined, features = c("CD19", "MS4A1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CD19", "MS4A1"))
# Cluster 14 (15 for them)

####RBCs & Platelets
VlnPlot(combined, features = c("HBA1", "PF4", "PPBP"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("PF4", "PPBP"))
# Cluster 15 (16 for them)

####Neutrophils
VlnPlot(combined, features = c("FCER1A", "CRLF2"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("FCER1A", "CRLF2"))
# Cluster 8 (8 for them()

####baso/eosinophils
VlnPlot(combined, features = c("CCR3", "PTGDR2"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CCR3", "PTGDR2"))
# Just set as 9 (their 10) because it matches their diagram

####Macrophages
VlnPlot(combined, features = c("CD14", "CD68"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CD14", "CD68"))
# Clusters 3,4,8,13,15 (16 for CD68)

#####M1
VlnPlot(combined, features = c("ITGAX", "HLA-DRA", "HLA-DRB1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("ITGAX", "HLA-DRA", "HLA-DRB1"))
# 3. Potentially 8 (3 and 14 for them)
# 3,4,8,13,15

#####M2
VlnPlot(combined, features = c("CD163", "MSR1", "MRC1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("CD163", "MSR1", "MRC1"))
# 3 and 8 and maybe 15 (5 and 13 for them)

# Endothelial cells
VlnPlot(combined, features = c("VWF", "TIE1"), slot = "counts", log = TRUE)
FeaturePlot(combined, features = c("VWF", "TIE1"))
# None, consistent with the paper

###########################################################
##Re-cluster by localization of cell-specific markers
new.cluster.ids <- list(
  "0" = "GC",
  "1" = "T-helper",
  "2" = "GC",
  "3" = "M2-Macrophage",
  "4" = "M1-Macrophage",
  "5" = "NKT",
  "6" = "NK",
  "7" = "Cytotoxic_T",
  "8" = "Neutrophil",
  "9" = "Baso_eosinophil",
  "10" = "GC",
  "11" = "GC",
  "12" = "GC",
  "13" = "M1-Macrophage",
  "14" = "B",
  "15" = "RBC_Platelet",
  "16" = "Dendritic",
  "17" = "TC"
)

# Add cell type names to metadata
combined <- RenameIdents(combined, new.cluster.ids)
combined[["cell_type"]] <- Idents(combined)

# Plot
DimPlot(combined, reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 5) + NoLegend()

# Remove scale data (there is an issue with seurat disk and the raw counts when converting to h5ad if not) - https://github.com/mojaveazure/seurat-disk/issues/75
seurat_ob <- DietSeurat(combined)

# https://github.com/mojaveazure/seurat-disk/issues/23
i <- sapply(seurat_ob@meta.data, is.factor)
seurat_ob@meta.data[i] <- lapply(seurat_ob@meta.data[i], as.character)

# Save as a h5ad file
SaveH5Seurat(seurat_ob, filename = out_seurat)
Convert(out_seurat, dest = "h5ad")


```
- Cells showing a minimum expression of 500 genes were filtered, and 
- genes showing the expression in a minimum of five cells were applied.  # it says 3 in the code
- After log-normalization using NormalizeData(), 
- four datasets were integrated following FindIntegrationAnchors() and IntegrateData() 
- procedures with the first 20 principal components for weighting. Filtered cells were dimensionally reduced and visualized on a Uniform Manifold Approximation and Projection (UMAP).25 Unsupervised clustering analyses on the UMAP embedding were conducted using #FindClusters() with “resolution = 0.05,” and cells were divided by four individual IVF patient samples in the same UMAP.

Women (age: 30–38 years) exhibiting regular menstrual cycles and had not taken hormonal contraceptives for at least 3 months before their enrollment in the study underwent laparoscopic sterilization.


- Follicular cell			CYP11A1, HSD3B2
- Granulosa cell (GC)		CYP19A1, INHA
- Theca cell (TC)	COL1A1, COL1A2, COL3A1,
- Endothelial cell		 VWF, TIE1
- Leukocytes			PTPRC
- Macrophage		CD14, CD68
- M1-macrophage		ITGAX, HLA-DRA, HLA-DRB1
- M2-macrophage		CD163, MSR1, MRC1
- Helper/cytotoxic T cell	CD3E, CD8A, CD4
- NK cell			NKG7, NCR1
- Neutrophil			FCER1A, CRLF2 
- Baso/eosinophil		CCR3
- B cell			CD19, MS4A1
- Dendritic cell		FCER1A 
- RBC			HBA1 
- Platelet	 	 	PF4, PPBP
