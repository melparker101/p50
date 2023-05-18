We have separate files for each sample, e.g. "GSM5761690_filtered_feature_bc_matrix_Don1.tar.gz". Unzip and extract from tar. They also provide R code on GEO.
Supplementary file: https://static-content.springer.com/esm/art%3A10.1038%2Fs42003-022-04384-8/MediaObjects/42003_2022_4384_MOESM2_ESM.pdf 
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

# "After routine quality control1, 48,147 cells remained for downstream for analysis."
# We only have 5607. Try reading in separately

################################################
# Load libraries
library(Seurat)

# Define data path
path = "data/counts/GSE192722"

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

VlnPlot(merged_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol = 4)

# QC
merged_ob[["percent.mito"]] <- PercentageFeatureSet(merged_ob, pattern = "^MT-")
merged_ob[["percent.ribo"]] <- PercentageFeatureSet(merged_ob, pattern = "^RP[SL]")

VlnPlot(merged_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), pt.size = 0, ncol = 4)

merged_ob <- subset(merged_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 20)

# 48147

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

# Histogram of library contribution to each cluster
# Idents(Theca) <- "seurat_clusters"
Theca[["seurat_clusters"]] <- Idents(Theca)
Clus0 <- subset(Theca, idents = c("0"))
Clus1 <- subset(Theca, idents = c("1"))
Clus2 <- subset(Theca, idents = c("2"))
Clus3 <- subset(Theca, idents = c("3"))
Clus4 <- subset(Theca, idents = c("4"))
Clus5 <- subset(Theca, idents = c("5"))
Clus6 <- subset(Theca, idents = c("6"))
Clus7 <- subset(Theca, idents = c("7"))
Clus8 <- subset(Theca, idents = c("8"))
Clus9 <- subset(Theca, idents = c("9"))
Clus10 <- subset(Theca, idents = c("10"))
Clus11 <- subset(Theca, idents = c("11"))
Clus12 <- subset(Theca, idents = c("12"))
Clus13 <- subset(Theca, idents = c("13"))
Clus14 <- subset(Theca, idents = c("14"))
Clus15 <- subset(Theca, idents = c("15"))
Clus16 <- subset(Theca, idents = c("16"))
Clus17 <- subset(Theca, idents = c("17"))
Clus18 <- subset(Theca, idents = c("18"))
Clus19 <- subset(Theca, idents = c("19"))
Clus20 <- subset(Theca, idents = c("20"))
Clus21 <- subset(Theca, idents = c("21"))
table(Clus0@meta.data$labels)
table(Clus1@meta.data$labels)
table(Clus2@meta.data$labels)
table(Clus3@meta.data$labels)
table(Clus4@meta.data$labels)
table(Clus5@meta.data$labels)
table(Clus6@meta.data$labels)
table(Clus7@meta.data$labels)
table(Clus8@meta.data$labels)
table(Clus9@meta.data$labels)
table(Clus10@meta.data$labels)
table(Clus11@meta.data$labels)
table(Clus12@meta.data$labels)
table(Clus13@meta.data$labels)
table(Clus14@meta.data$labels)
table(Clus15@meta.data$labels)
table(Clus16@meta.data$labels)
table(Clus17@meta.data$labels)
table(Clus18@meta.data$labels)
table(Clus19@meta.data$labels)
table(Clus20@meta.data$labels)
table(Clus21@meta.data$labels)

# Subset Theca Cells
ThecaStroma <- subset(Theca, idents = c("1","2","3","4","5","12"))
# 24033 cells. They have 23,736 cells for this






```
