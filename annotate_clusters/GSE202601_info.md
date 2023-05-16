Paper: https://www.biorxiv.org/content/biorxiv/early/2022/05/19/2022.05.18.492547.full.pdf


```
library(data.table)
library(Seurat)
library(magrittr)
library(cowplot)
library(harmony)
library(SeuratData)

# Dataset
dataset <- "GSE202601"

# Data path
path <- paste("data/counts",dataset,sep="/")

# Input files
mat_file <- paste(path,"GSE202601_human_ovary_snRNA-seq_count.rds",sep="/")
meta_file <- paste(path,"GSE202601_human_ovary_snRNA-seq_metadata.txt",sep="/")

# Load in data
mat <- readRDS(file = mat_file)
meta <- as.data.frame(fread(meta_file, header=T))

colnames(meta)[1] <- "cell_id"

# See 'Metadata format' in https://github.com/perslab/CELLEX
meta <- meta[,c(1,5)]
colnames(meta) <- c("cell_id","abbr")

###
# Format metadata
rownames(meta) <- meta$cell_id
meta <- meta[,-1]

# Create Seurat object
ovary <- CreateSeuratObject(counts = mat, project = "ovary", min.cells = 3, min.features = 200, meta.data = meta)

#######################################
# Code chunk from original author
ovary[["percent.mt"]] <- PercentageFeatureSet(ovary, pattern = "^MT-")
ovary<- subset(ovary, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
ovary<- NormalizeData(ovary, normalization.method = "LogNormalize", scale.factor = 10000)
ovary<- FindVariableFeatures(ovary, selection.method = "vst", nfeatures = 2100)
all.genes <- rownames(ovary)
ovary<- ScaleData(ovary, features = all.genes)
ovary<- RunPCA(ovary, features = VariableFeatures(object = ovary))

######################################

# Code from original authors
library(Seurat)
library(dplyr)
library(harmony)
ovary_count<-read.table("ovary_count",header=T,row.name=1,sep="\t")
ovary<- CreateSeuratObject(counts = ovary_count, project = "ovary")

# metadata <- meta

ovary[["percent.mt"]] <- PercentageFeatureSet(ovary, pattern = "^MT-")

ovary<- NormalizeData(ovary, normalization.method = "LogNormalize", scale.factor = 10000)
ovary<- FindVariableFeatures(ovary, selection.method = "vst", nfeatures = 2100)
all.genes <- rownames(ovary)
ovary<- ScaleData(ovary, features = all.genes,vars.to.regress = "percent.mt")
ovary<- RunPCA(ovary, features = VariableFeatures(object = ovary))

ovary@metadata$age<-ovary_metadata$age
ovary@metadata$group<-ovary_metadata$group
ovary <- ovary %>% RunHarmony("age", plot_convergence = TRUE)
ovary <- ovary %>% RunUMAP(reduction = "harmony", dims = 1:15) %>% FindNeighbors(reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.1) %>% identity()

DimPlot(ovary,pt.size = .5)
DimPlot(ovary,pt.size = .5,split.by="group")
DimPlot(ovary,pt.size = .5,split.by="age")

ovary_roc.markers <- FindAllMarkers(ovary, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "roc")

new.cluster.ids <- c("Stromal_cell", "Blood_endothelial_cell", "Granulosa_cell", "Smooth_muscle_cell", "Immune_cell", "Lymphatic_endothelial_cell", "Epithelial_cell", "Theca_cell","Stromal_cell")
names(new.cluster.ids) <- levels(ovary)

ovary <- RenameIdents(ovary, new.cluster.ids)
DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5,split.by="age")
DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.5,split.by="group")

DotPlot(ovary,features = c("DCN","COL6A3","LUM","PDGFRA","VWF","FLT1","CDH2","SERPINE2","CYP19A1","INHA","FOXL2","AMH","MYH11","ACTA2","TAGLN","MCAM","PTPRC","CD53","CXCR4","PROX1","FLT4","PAX8","CDH1","CLDN1","STAR","CYP17A1"),cols = c("white","red"))

SC_aging_markers <- FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Stromal_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
BEC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Blood_endothelial_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
GC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Granulosa_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
SMC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Smooth_muscle_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
IC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Immune_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
LEC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Lymphatic_endothelial_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
EpiC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Epithelial_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")
TC_aging_markers<-FindMarkers(ovary, ident.1 = "old", ident.2 ="young",subset.ident = "Theca_cell",group.by="group", min.pct = 0.25,logfc.threshold=0.25,test.use="MAST")


SC_aging_logfc<-FindMarkers(ovary,subset.ident = "Stromal_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
BEC_aging_logfc<-FindMarkers(ovary,subset.ident = "Blood_endothelial_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
GC_aging_logfc<-FindMarkers(ovary,subset.ident = "Granulosa_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
SMC_aging_logfc<-FindMarkers(ovary,subset.ident = "Smooth_muscle_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
IC_aging_logfc<-FindMarkers(ovary,subset.ident = "Immune_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
LEC_aging_logfc<-FindMarkers(ovary,subset.ident = "Lymphatic_endothelial_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
EpiC_aging_logfc<-FindMarkers(ovary,subset.ident = "Epithelial_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")
TC_aging_logfc<-FindMarkers(ovary,subset.ident = "Theca_cell",ident.1 = "old",ident.2 = "young",features = all.genes,group.by = "group",logfc.threshold = -Inf,min.pct = -Inf,test.use = "MAST")

SC_aging_logfc<-SC_aging_logfc[,"avg_log2FC",drop=F]
BEC_aging_logfc<-BEC_aging_logfc[,"avg_log2FC",drop=F]
GC_aging_logfc<-GC_aging_logfc[,"avg_log2FC",drop=F]
SMC_aging_logfc<-SMC_aging_logfc[,"avg_log2FC",drop=F]
IC_aging_logfc<-IC_aging_logfc[,"avg_log2FC",drop=F]
LEC_aging_logfc<-LEC_aging_logfc[,"avg_log2FC",drop=F]
EpiC_aging_logfc<-EpiC_aging_logfc[,"avg_log2FC",drop=F]
TC_aging_logfc<-TC_aging_logfc[,"avg_log2FC",drop=F]
ovary_aging_logfc<-merge(SC_aging_logfc,BEC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,GC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,SMC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,IC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,LEC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,EpiC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
ovary_aging_logfc<-merge(ovary_aging_logfc,TC_aging_logfc,by=0,all=F)
rownames(ovary_aging_logfc)<-ovary_aging_logfc$Row.names
ovary_aging_logfc<-ovary_aging_logfc[,-1]
colnames(ovary_aging_logfc)<-c("SC","BEC","GC","SMC","IC","LEC","EpiC","TC")

#######################################

# Plot 
DimPlot(ovary, group.by = c("new.cluster.ids"), ncol = 2, label = TRUE)




```
