![QC](https://github.com/melparker101/p50_Infertility/assets/98864236/58c6859c-6eb1-41d2-9d1b-3635806cba10)

QC:
- 200–5000 genes, UMIs less than 30,000, and a mitochondrial gene expression percentage of less than 50% for the following analyses.

```
####################################

library(data.table)
library(Seurat)
library(magrittr)
library(cowplot)
library(harmony)
library(SeuratData)
library(SeuratDisk)
library(ggplot2)

dataset <- "GSE189960"
path <- "data/counts/GSE189960"

out_seurat <- paste0(path, "/", dataset, "_counts.h5Seurat")

filelist <- list.files(path,pattern='.expression_matrix.txt')
patient_id <- c()

for(file in filelist){
  # Read in data and create seurat_ob
  dt <- fread(paste(path,file,sep="/"))
  df <- as.data.frame(dt)
  rownames(df) <- df$V1
  df <- df[,-1]
  name <- gsub(".*?([[:alnum:]]+)[.].*", "\\1", file)
  seurat_ob <- CreateSeuratObject(counts = df, project=name, min.cells = 3, min.features = 200)
  
  # Rename Seurat object
  assign(name, seurat_ob)
  rm(seurat_ob)
  patient_id <- append(patient_id,name)
}

# 14,592

# print(max(seurat_ob@meta.data$nCount_RNA))

# Merge Seurat objects. The line below is what we want to achieve, but we make it more reproducable
# merged_ob <- merge(P1, y= c(P2,P3,P4,P5,P6), add.cell.ids = patient_id)

# Create a command for merging seurat objects
# Create a list of P2-P5 separated by commas
merge_list <- paste(patient_id[2:length(patient_id)], collapse = ",")

# "add.cell.ids = patient_id" to avoid duplicate cell names
mergeSamples <- paste0("merged_ob <- merge(",patient_id[1],", y = c(", merge_list, "), add.cell.ids = patient_id, project = \"ovaries,min.cells=3 \")")

# Run command to merge Seurat objects
merged_ob <- eval(parse(text = mergeSamples))

merged_ob[["percent.mt"]] <- PercentageFeatureSet(merged_ob, pattern = "^MT-")
VlnPlot(merged_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "orig.ident")

merged_ob <- subset(merged_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 50 & nCount_RNA < 30000)

table(Idents(merged_ob))

# Add a sample column
merged_ob[["sample"]] <- Idents(object = merged_ob)

# This plot will be different to the original plot in supplementary figure 1 because theirs is produced before filtering
VlnPlot(merged_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "orig.ident")

# We retained cells with 200–5000 genes, UMIs less than 30,000, and a mitochondrial gene expression percentage of less than 50% for the following analyses. 
# We end up with 14,857 cells rather than 14,592 cells. I am not sure why and they do not provide extra filtering information, but it should not make a significant difference.

# NormalizeData and ScaleData to normalize and scale the gene expression matrix. For the PCA analysis, the top 2000 variable genes were chosen by FindVariableFeatures

# Clone the Seurat object to see how it clusters without batch effect removal
# merged_ob_noBER <- merged_ob

# Without batch effect removal
# merged_ob_noBER <- merged_ob_noBER %>%
#     NormalizeData() %>%
#     FindVariableFeatures(nfeatures = 2000) %>%
#     ScaleData() %>%
#     RunPCA() %>%
#     RunUMAP(reduction = "pca", dims = 1:20) %>%
#     FindNeighbors(dims = 1:20) %>%
#     FindClusters(resolution = 1.2)
#     # RunTSNE()
#     # DimPlot(reduction = "tsne")

# Visualise the UMAP
# p1 <- DimPlot(merged_ob_noBER, reduction = "umap", group.by = "sample")
# p2 <- DimPlot(merged_ob_noBER, reduction = "umap", label = TRUE)
# plot_grid(p1, p2)
# The batch effect is massive

# With batch effect removal using harmony
merged_ob <- merged_ob %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = "sample") %>%
    RunUMAP(reduction = "harmony", dims = 1:20) %>%
    # RunUMAP(reduction = "pca", dims = 1:20) %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    # FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1.2)
    # RunTSNE()
    # DimPlot(reduction = "tsne")

# Add the clusters as a column
table(Idents(object = merged_ob))
merged_ob[["seurat_clusters"]] <- Idents(object = merged_ob)

p1 <- DimPlot(merged_ob, reduction = "umap", group.by = "sample")
p2 <- DimPlot(merged_ob, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
plot_grid(p1, p2)

# This does the same
# DimPlot(merged_ob, group.by = c("sample", "seurat_clusters"), ncol = 2, label = TRUE)

DimPlot(merged_ob, group.by = "seurat_clusters", label = TRUE)

FeaturePlot(merged_ob, features = c("STAR", "CD68", "CD1C", "FCGR3B", "CD3D", "CDH1", "SERPINE2", "CD163",
    "FCER1A", "CXCR2", "CD3G", "EPCAM"))
    
FeaturePlot(merged_ob, features = "SERPINE1")

# DimPlot(object = merged_ob, reduction = "pca")
# DimPlot(object = merged_ob, reduction = "umap")

# Cells were separated into 23 clusters by FindClusters, by using the top 20 principle components and a resolution parameter of 1.2. 
# For the clustering of GCs and macrophages, we set the resolution to 1.2 and applied the uniform mainfold approximation and projection (UMAP) algorithm to visualize cells in a two-dimensional space

# We can manually annotate... 
# 9 GC clusters: 0,2,3,4,5,8,9,10,14
# 5 macrophages?:1,6,7,11,13,15,16,19,21 - figure this out
# DCs:
# T cells:
# neutrophils:

FeaturePlot(merged_ob, features = c("STAR", "SERPINE2"))  # Granulosa
FeaturePlot(merged_ob, features = c("PTPRC", "CD68", "CD163")) # Macrophage
FeaturePlot(merged_ob, features = c("PTPRC", "CD1C", "FCER1A"))  # DC
FeaturePlot(merged_ob, features = c("PTPRC", "CXCR2", "FCGR3B"))  # Neutrophil
FeaturePlot(merged_ob, features = c("PTPRC", "CD3D", "CD3E", "CD3G"))  # T_cell
FeaturePlot(merged_ob, features = c("EPCAM", "KRT18", "CDH1"))  # Epithelium

# Cell labels
# GC_clusters <- c("0","1","2","5","6","9","11","14","15","18","21")
# neutrophil_clusters <- "16"
# macrophage_clusters <- c("3","4","8","13","17","19")
# Tcell_clusters <- c("12","22")
# epithelium_clusters <- c("10","20")
# DC_clusters <- "7"

# Make a dictionary for the cell labels
cluster_dict <- c(
  "0" = "GC",
  "1" = "GC",
  "2" = "GC",
  "5" = "GC",
  "6" = "GC",
  "9" = "GC",
  "11" = "GC",
  "14" = "GC",
  "15" = "GC",
  "18" = "GC",
  "21" = "GC",
  "16" = "neutrophil",
  "3" = "macrophage",
  "4" = "macrophage",
  "8" = "macrophage",
  "13" = "macrophage",
  "17" = "macrophage",
  "19" = "macrophage",
  "12" = "Tcell",
  "22" = "Tcell",
  "10" = "epithelium",
  "20" = "epithelium",
  "7" = "DC"
)

# Reorder dictionary based on clusters
cluster_dict <- cluster_dict[as.character(sort(as.numeric(names(cluster_dict))))]

# Rename the idents and add as a column in seurat object
merged_ob <- RenameIdents(merged_ob, cluster_dict)
merged_ob[["cell_type"]] <- Idents(merged_ob)

# Flip x and y axis and plot umap
DimPlot(merged_ob, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() + scale_x_reverse() + scale_y_reverse()

# Remove scale data (there is an issue with seurat disk and the raw counts when converting to h5ad if not) - https://github.com/mojaveazure/seurat-disk/issues/75
seurat_ob <- DietSeurat(merged_ob)

# https://github.com/mojaveazure/seurat-disk/issues/23
i <- sapply(seurat_ob@meta.data, is.factor)
seurat_ob@meta.data[i] <- lapply(seurat_ob@meta.data[i], as.character)

# Save as a h5ad file
SaveH5Seurat(seurat_ob, filename = out_seurat)
Convert(out_seurat, dest = "h5ad")

```
Canonical markers and highly differentially expressed genes (DEGs) enabled us to identify six major cell types: 
- STAR+SERPINE2+ granulosa cells (GCs, 9614 cells, 66%) [44,45], 
- PTPRC+CD68+ macrophages (3601 cells, 25%), 
- PTPRC+CD1C+ dendritic cells (DCs, 754 cells, 5%), 
- PTPRC+CXCR2+ neutrophils (168 cells, 1%), 
- PTPRC+CD3D+CD3E+ T cells (260 cell, 2%)
- EPCAM+KRT18+ epithelial cells (195 cells, 1%) (Figure 1D–F)

Using seurat v 4.3 and manually picking based on the above marker genes and comparing to the original cluster graph from the authors:
16 - Neutrophil
0,1,2,5,6,9,11,14,15,18,21 - GC
3,4,8,13,17,19 - Macrophages
7 - DC
12,22 - T cells
10,20 - Epithelium





