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

path <- "data/counts/GSE189960"

filelist <- list.files(path,pattern='.expression_matrix.txt')
patient_id <- c()

for(file in filelist){
  # Read in data and create seurat_ob
  dt <- fread(paste(path,file,sep="/"))
  df <- as.data.frame(dt)
  rownames(df) <- df$V1
  df <- df[,-1]
  name <- gsub(".*?([[:alnum:]]+)[.].*", "\\1", file)
  seurat_ob <- CreateSeuratObject(counts = df, project=name)
  
  # QC
  seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")
  seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 50 & nCount_RNA < 30000)
  
  # Rename Seurat object
  assign(name, seurat_ob)
  rm(seurat_ob)
  patient_id <- append(patient_id,name)
}

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

table(Idents(merged_ob))

# Add a sample column
merged_ob[["Sample"]] <- Idents(object = merged_ob)

# This plot will be different to the original plot in supplementary figure 1 because theirs is produced before filtering
VlnPlot(merged_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "orig.ident")

# We retained cells with 200–5000 genes, UMIs less than 30,000, and a mitochondrial gene expression percentage of less than 50% for the following analyses. 
# We end up with 14,857 cells rather than 14,592 cells. I am not sure why and they do not provide extra filtering information, but it should not make a significant difference.

NormalizeData
and ScaleData to normalize and scale the gene expression matrix. For the PCA analysis, the top 2000 variable genes were chosen by FindVariableFeatures

# merged_ob2 <- merged_ob
merged_ob <- merged_ob2

<!-- 
merged_ob <- NormalizeData(object = merged_ob)
merged_ob <- FindVariableFeatures(object = merged_ob, nfeatures = 2000)
merged_ob <- ScaleData(object = merged_ob)
merged_ob <- RunPCA(object = merged_ob)
merged_ob <- RunUMAP(object = merged_ob, reduction = "pca", dims= 1:20)
merged_ob <- FindNeighbors(object = merged_ob,dims = 1:20)
merged_ob <- FindClusters(object = merged_ob,resolution = 1.2) -->


merged_ob_h <- merged_ob2
merged_ob_h[["Sample"]] <- Idents(object = merged_ob_h)

# Without batch effect removal
merged_ob <- merged_ob %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunUMAP(reduction = "pca", dims = 1:20) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1.2)
    # RunTSNE()
    # DimPlot(reduction = "tsne")

# Visualise the UMAP
p1 <- DimPlot(merged_ob, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(merged_ob, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
# The batch effect is massive

# With batch effect removal using harmony
merged_ob_h <- merged_ob_h %>%
    NormalizeData() %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = "Sample") %>%
    RunUMAP(reduction = "harmony", dims = 1:20) %>%
    # RunUMAP(reduction = "pca", dims = 1:20) %>%
    FindNeighbors(reduction = "harmony", dims = 1:20) %>%
    # FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 1.2)
    # RunTSNE()
    # DimPlot(reduction = "tsne")

table(Idents(object = merged_ob_h))
merged_ob_h[["clusters"]] <- Idents(object = merged_ob_h)

p1 <- DimPlot(merged_ob_h, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(merged_ob_h, reduction = "umap", group.by = "clusters", label = TRUE)
plot_grid(p1, p2)

DimPlot(merged_ob_h, group.by = c("Sample", "clusters"), ncol = 2, label = TRUE)

FeaturePlot(merged_ob_h, features = c("STAR", "CD68", "CD1C", "FCGR3B", "CD3D", "CDH1", "SERPINE2", "CD163",
    "FCER1A", "CXCR2", "CD3G", "EPCAM"))
    
 FeaturePlot(merged_ob_h, features = "SERPINE1")

# DimPlot(object = merged_ob, reduction = "pca")
# DimPlot(object = merged_ob, reduction = "umap")

# Cells were separated into 23 clusters by FindClusters, by using the top 20 principle components and a resolution parameter of 1.2. 
# For the clustering of GCs and macrophages, we set the resolution to 1.2 and applied the uniform mainfold approximation and projection (UMAP) algorithm to visualize cells in a two-dimensional space

# DimPlot(object = merged_ob, reduction = "pca")
# DimPlot(object = merged_ob, reduction = "umap")

# We can manually annotate... 
# 9 GC clusters: 0,2,3,4,5,8,9,10,14
# 5 macrophages?:1,6,7,11,13,15,16,19,21 - figure this out
# DCs:
# T cells:
# neutrophils:

FeaturePlot(merged_ob_h, features = c("STAR", "SERPINE2"))
FeaturePlot(merged_ob_h, features = c("PTPRC", "CD68"))

```
Canonical markers and highly differentially expressed genes (DEGs) enabled us to identify six major cell types: 
- STAR+SERPINE2+ granulosa cells (GCs, 9614 cells, 66%) [44,45], 
- PTPRC+CD68+ macrophages (3601 cells, 25%), 
- PTPRC+CD1C+ dendritic cells (DCs, 754 cells, 5%), 
- PTPRC+CXCR2+ neutrophils (168 cells, 1%), 
- PTPRC+CD3D+CD3E+ T cells (260 cell, 2%)
- EPCAM+KRT18+ epithelial cells (195 cells, 1%) (Figure 1D–F)

