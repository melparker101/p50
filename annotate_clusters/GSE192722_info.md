We have separate files for each sample, e.g. "GSM5761690_filtered_feature_bc_matrix_Don1.tar.gz". Unzip and extract from tar. They also provide R code on GEO.
Supplementary file: https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-022-04384-8/MediaObjects/42003_2022_4384_MOESM2_ESM.pdf 
``` R
# Load libraries
library(Seurat)
library(SeuratDisk)
library(ggplot2)

# Define dataset, data path and seurat outfile
dataset <- "GSE192722"
path <- paste0("data/counts/",dataset)
out_seurat <- paste0(path, "/", dataset, "_counts.h5Seurat")

# Make a list of the data directories for the 8 different samples
data_dirs <- list.files(path, "(Don|Pt)")
data_dirs <- data_dirs[!grepl(".tar", data_dirs)]

# Make an empty vector to add sample names to
samples <- c()

# Loop through samples
for(dir in paste(path,data_dirs, sep="/")){
  # Extract sample name
  name <- gsub("data/counts/GSE192722/filtered_feature_bc_matrix_","",dir)
  # Read in data
  data <- Read10X(data.dir = dir)
  # Create Seurat object
  seurat_object <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 200)
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

# Now we have 63301 features. Filter.
# 62255

# >200 but <4000 features and <20% of reads mapping to mitochondrial genes were retained.
# They actually filtered for <7500 features
# table(Theca@meta.data$percent.mito < 0.20 & Theca@meta.data$nFeature_RNA > 200 & Theca@meta.data$nFeature_RNA < 7500)

# QC
# Add mitochondria percent and ribosomal percent columns
merged_ob[["percent.mito"]] <- PercentageFeatureSet(merged_ob, pattern = "^MT-")
merged_ob[["percent.ribo"]] <- PercentageFeatureSet(merged_ob, pattern = "^RP[SL]")

# Plot
VlnPlot(merged_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), pt.size = 0, ncol = 4)

# Filter
merged_ob <- subset(merged_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 20)
# 48147

# Plot again
VlnPlot(merged_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), pt.size = 0, ncol = 4)

# Checkpoint
Theca <- merged_ob

# Batch effect?
# Idents(Theca) <- "labels"
Don1_5 <- WhichCells(Theca, idents = c("Don1", "Don2", "Don3", "Don4", "Don5"))
Pt2 <- WhichCells(Theca, idents = c("Pt2"))
Don6_11 <- WhichCells(Theca, idents = c("Don6_11", "Pt1"))

Theca[["samples"]] <- Idents(Theca)
Theca <- SetIdent(Theca, Don1_5, "Don1_5")
Theca <- SetIdent(Theca, Pt2, "Pt2")
Theca <- SetIdent(Theca, Don6_11, "Don6_11")
table(Idents(Theca))
# Theca@meta.data$batch <- Theca@active.ident
# Idents(Theca) <- "batch"
Theca[["batch"]] <- Idents(Theca)

# Split the dataset into a list of two seurat objects () # It should be 3?
Theca.list <- SplitObject(Theca, split.by = "batch")

# Normalise and identify variable features for each dataset independently
Theca.list <- lapply(X = Theca.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Theca.list)
Theca.list <- lapply(X = Theca.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
Theca.anchors <- FindIntegrationAnchors(object.list = Theca.list, anchor.features = 2000, reduction = "rpca", dims = 1:50)
Theca <- IntegrateData(anchorset = Theca.anchors, dims = 1:50)
Theca <- ScaleData(Theca, verbose = FALSE)

# Plot According to Cell Cycle
Theca <- CellCycleScoring(Theca, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
# Regress Cell Cycle Effects (if desired)
Theca <- ScaleData(Theca, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Theca))

Theca <- RunPCA(object = Theca)
Theca <- FindNeighbors(object = Theca)
Theca <- FindClusters(object = Theca)
Theca <- RunUMAP(Theca, dims = 1:50)

# Plots
# DimPlot(Theca, reduction = "umap", group.by = "labels", pt.size = .05, label = TRUE)
DimPlot(Theca, reduction = "umap", group.by = "samples", pt.size = .05, label = TRUE)
DimPlot(Theca, reduction = "umap", group.by = "samples", pt.size = .05, label = TRUE)+ NoLegend()
DimPlot(Theca, reduction = "umap", group.by = "seurat_clusters", pt.size = .05, label = TRUE)+ NoLegend()
DimPlot(Theca, reduction = "umap", group.by = "Phase", pt.size = .05, label = TRUE)+ NoLegend()

DimPlot(Theca, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

FeaturePlot(object = Theca, features = "PTPRC", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)
FeaturePlot(object = Theca, features = "CDH5", reduction = "umap",cols = c("white","red"),order = TRUE,pt.size = 0.25)
FeaturePlot(object = Theca, features = "RGS5", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)
FeaturePlot(object = Theca, features = "AMH", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)
FeaturePlot(object = Theca, features = "CYP17A1", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)
FeaturePlot(object = Theca, features = "OGN", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)
FeaturePlot(object = Theca, features = "ACTA2", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)
FeaturePlot(object = Theca, features = "KRT8", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)
FeaturePlot(object = Theca, features = "MCAM", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)

# Plot for theca cells
FeaturePlot(object = Theca, features = "AMH", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)  # GC
FeaturePlot(object = Theca, features = "PTPRC", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)  # HEM
FeaturePlot(object = Theca, features = "CDH5", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)  # EC
FeaturePlot(object = Theca, features = "RGS5", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)  # SMC
FeaturePlot(object = Theca, features = "TCF21", reduction = "umap",cols = c("grey90","red"),order = TRUE,pt.size = 0.25)  # Theca/stroma ?
# OE? KRT8? PTGDS?

# VlnPlot(object = ThecaProp, features = "TCF21", pt.size = 0, group.by = "seurat_clusters") +
#   stat_summary(fun = median, geom='point', size = 25, colour = "black", shape = 95)+NoLegend()

# Subset Theca Cells
ThecaStroma <- subset(Theca, idents = c("1","2","3","4","5","12"))
# 24033 cells. They have 23,736 cells for this

# Assign cluster numbers to a cell type
Theca_stroma_clusters <- c("1","2","3","4","5","12")
OE_clusters <- c("6","11","16")
GC_clusters <- c("7","15")
SMC_clusters <- c("14","20")
EC_clusters <- c("9","13")
HEM_clusters <- c("17","18")

# Create the cluster_dict dictionary
cluster_dict <- c(
  setNames(rep("Theca_stroma", length(Theca_stroma_clusters)), Theca_stroma_clusters),
  setNames(rep("OE", length(OE_clusters)), OE_clusters),
  setNames(rep("GC", length(GC_clusters)), GC_clusters),
  setNames(rep("SMC", length(SMC_clusters)), SMC_clusters),
  setNames(rep("EC", length(EC_clusters)), EC_clusters),
  setNames(rep("HEM", length(HEM_clusters)), HEM_clusters)
)

# Reorder dictionary based on clusters
cluster_dict <- cluster_dict[as.character(sort(as.numeric(names(cluster_dict))))]

# Rename the idents and add as a column in seurat object
Theca <- RenameIdents(Theca, cluster_dict)
Theca[["cell_type"]] <- Idents(Theca)

# Flip and plot umap
DimPlot(Theca, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + scale_x_reverse()

# Make a list of cluster names to keep
clusters_keep <- unique(as.vector(cluster_dict))

# Subset to only include the labelled clusters. Rename the object again so we can reuse code from another dataset when saving
merged_ob <- subset(x = Theca, idents = clusters_keep)

# Remove scale data (there is an issue with seurat disk and the raw counts when converting to h5ad if not) - https://github.com/mojaveazure/seurat-disk/issues/75
seurat_ob <- DietSeurat(merged_ob)

# https://github.com/mojaveazure/seurat-disk/issues/23
i <- sapply(seurat_ob@meta.data, is.factor)
seurat_ob@meta.data[i] <- lapply(seurat_ob@meta.data[i], as.character)

# Save as a h5ad file
SaveH5Seurat(seurat_ob, filename = out_seurat)
Convert(out_seurat, dest = "h5ad")

# Theca = Seurat object with all clusters
# merged_ob = Seurat object with only annotated cell types in
# seurat_ob = Seurat object with scale data removed and metadata factor type changed to character
```
