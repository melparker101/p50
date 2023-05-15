![QC](https://github.com/melparker101/p50_Infertility/assets/98864236/58c6859c-6eb1-41d2-9d1b-3635806cba10)

QC:
- 200–5000 genes, UMIs less than 30,000, and a mitochondrial gene expression percentage of less than 50% for the following analyses.

```
library(data.table)
library(Seurat)

# Read in counts and create a Seurat object
dt <- fread(paste(path,"GSM5710585_P1.expression_matrix.txt",sep="/"))
df <- as.data.frame(dt)
rownames(df) <- df$V1
df <- df[,-1]
seurat_ob <- CreateSeuratObject(counts = df)
seurat_ob  # 1786 cells

# Follow their QC
# seurat_ob <- CreateSeuratObject(counts = df, project = "ovaries", min.cells = 3, min.features = 200)
seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")

# Visualise
VlnPlot(seurat_ob, features=c("nCount_RNA","percent.MT", "percent.Ribosomal","percent.Largest.Gene"))
VlnPlot(seurat_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(seurat_ob, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter
seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 50)

# For non-UMI data, nUMI represents the sum of the non-normalized values within a cell 
# nCount_RNA is UMIs/cell.

max(seurat_ob@meta.data$nCount_RNA)  # 27,160

####################################

library(data.table)
library(Seurat)

path <- "data/counts/GSE189960"

filelist <- list.files(path,pattern='.expression_matrix.txt')

for(file in filelist){
  # Read in data and create seurat_ob
dt <- fread(paste(path,file,sep="/"))
  df <- as.data.frame(dt)
  rownames(df) <- df$V1
  df <- df[,-1]
  seurat_ob <- CreateSeuratObject(counts = df)
  
  # QC
  seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")
  seurat_ob <- subset(seurat_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 50 & nCount_RNA < 30000)
  print(max(seurat_ob@meta.data$nCount_RNA))  # 27,160
    
  name <- gsub(".*?([[:alnum:]]+)[.].*", "\\1", file)
  assign(name, seurat_ob)
  rm(seurat_ob)
}




```
