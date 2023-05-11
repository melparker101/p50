#########################################################
# This code is for annotating clusters following the tutorial in the link below
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_06_celltype.html#scPred
# melodyjparker14@gmail.com - May 23
# Use dataset GSE118127 as a reference
# Just do basic processing for now
#########################################################

############################
# Load packages
############################
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
    library(scPred)
    library(SeuratDisk)
})

# Load the data
alldata <- readRDS("data/results/covid_qc_dr_int_cl.rds")
ctrl = alldata[, alldata$orig.ident == "ctrl_13"]

# set active assay to RNA and remove the CCA assay
ctrl@active.assay = "RNA"
ctrl[["CCA"]] = NULL
ctrl

############################
# Load in ref data
############################
ref_acc <- "GSE118127"
dataset1_acc <- "GSE202601"
dataset2_acc <- "GSE213216"

ref_file <- paste0("data/counts/",ref_dataset,"/local.rds")
# out_file <- paste0("annotating_clusters/", ref_dataset,)

# Load reference in Seurat object
reference <- readRDS(file = ref_file)
reference

# Have a look at the metadata
reference[[c('cluster_id','cell_description','cell_type')]][1:20,]
table(reference[['cluster_id']])

# Check dimensions
dim(GetAssayData(object = reference))  # 32922 genes

############################
# Load in other data
############################
dataset1_acc <- "GSE202601"
dataset2_acc <- "GSE213216"

dataset1_file <- paste0("data/counts/",dataset1_acc,"/counts_GSE202601.h5Seurat")
dataset2_file <- paste0("data/counts/",dataset2_acc,"/local.rds")

# Load reference in Seurat object
dataset1 <- LoadH5Seurat(dataset1_file)
reference

# Check dimensions
dim(GetAssayData(object = reference))

############################
# Run analysis pipeline for reference
############################

# Check that the data is normalised
GetAssayData(object = reference)[1:10,1:15]
GetAssayData(object = reference, slot = "counts")[1:10,1:15]
# The data is already normalised

# Check if variable gene selection has been run
head(HVFInfo(object = reference))
# Produces error - HGV needs running

# Check if data has been scaled
reference[["RNA"]]@scale.data[1:10,1:15]
# It hasn't been scaled

# Dimensionality reduction
# There is already a PCA saved
DimPlot(reference, reduction = "pca")
# We can rerun though

reference <- seurat_ob

# Run all analysis steps in one go using the magittr package
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
reference <- reference %>%
    # NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    RunUMAP(dims = 1:30)

# Plot UMAP for different hierarche of labels
DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE) + NoAxes()
DimPlot(reference, group.by = "cell_description", label = TRUE, repel = TRUE) + NoAxes()

# Top 10 HVGs
top10 <- head(VariableFeatures(reference), 10)
top10

# Plot HVGs
plot1 <- VariableFeaturePlot(reference)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Examine and visualize PCA results a few different ways
print(reference[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(reference, dims = 1:2, reduction = "pca")
DimHeatmap(reference, dims = 1, cells = 500, balanced = TRUE)

############################
# Run analysis pipeline for other datasets
############################

# Look at data
dataset1
GetAssayData(object = dataset1)[1:10,1:15]

dataset1 <- dataset1 %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    RunUMAP(dims = 1:30)
    
DimPlot(dataset1, group.by = "cell_type", label = TRUE, repel = TRUE) + NoAxes()


