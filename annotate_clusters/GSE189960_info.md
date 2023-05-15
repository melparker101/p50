![QC](https://github.com/melparker101/p50_Infertility/assets/98864236/58c6859c-6eb1-41d2-9d1b-3635806cba10)

QC:
- 200–5000 genes, UMIs less than 30,000, and a mitochondrial gene expression percentage of less than 50% for the following analyses.

```
####################################

library(data.table)
library(Seurat)
library(magrittr)

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

# This plot will be different to the original plot in supplementary figure 1 because theirs is produced before filtering
VlnPlot(merged_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, split.by = "orig.ident")

# We retained cells with 200–5000 genes, UMIs less than 30,000, and a mitochondrial gene expression percentage of less than 50% for the following analyses. 
# We end up with 14,857 cells rather than 14,592 cells. I am not sure why and they do not provide extra filtering information, but it should not make a significant difference.

NormalizeData
and ScaleData to normalize and scale the gene expression matrix. For the PCA analysis, the top 2000 variable genes were chosen by FindVariableFeatures

merged_ob2 <- merged_ob
merged_ob <- merged_ob2

merged_ob <- NormalizeData(object = merged_ob)
merged_ob <- FindVariableFeatures(object = merged_ob, nfeatures = 2000)
merged_ob <- ScaleData(object = merged_ob)
merged_ob <- RunPCA(object = merged_ob)
merged_ob <- RunUMAP(object = merged_ob, reduction = "pca", dims= 1:20)
merged_ob <- FindNeighbors(object = merged_ob,dims = 1:20)
merged_ob <- FindClusters(object = merged_ob,resolution = 1.2)


DimPlot(object = merged_ob, reduction = "pca")
DimPlot(object = merged_ob, reduction = "umap")

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

GetAssayData(object = merged_ob)[1:10,1:10]

# Cells were separated into 23 clusters by FindClusters, by using the top 20 principle components and a resolution parameter of 1.2. 
# For the clustering of GCs and macrophages, we set the resolution to 1.2 and applied the uniform mainfold approximation and projection (UMAP) algorithm to visualize cells in a two-dimensional space


library(harmony)
library(Seurat)
library(SeuratData)

```
