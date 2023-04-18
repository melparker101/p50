
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(dplyr)

dataset <- "GSE213216"

in_file <- "counts/GSE213216/counts_GSE213216.h5Seurat"
out_dir <- "cluster_markers/GSE213216"

# Create output directory
dir.create(out_dir)

# Read in Seurat object (this is already filtered for the samples we want to use)
seurat_ob <- LoadH5Seurat(in_file)
seurat_ob

# View metadata
seurat_ob[[c('seurat_clusters','active_cluster','DF.classifications')]][1:20,]
table(seurat_ob[['seurat_clusters']])
table(seurat_ob[['active_cluster']])

# Remove doublets
table(seurat_ob[['DF.classifications']])
doub_rem <- subset(seurat_ob, subset=DF.classifications %in% "Singlet")

##### Find Markers #####

### 1. Use active clusters (9 clusters)
clust_no <- "9C"

# Set indentity classes as the pre-calculated Seurat clusters
Idents(object = seurat_ob) <- "active_cluster"
levels(seurat_ob)

# Find markers for every cluster compared to all remaining cells, report only the positive ones
combined_markers <- FindAllMarkers(object = seurat_ob, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)  
View(combined_markers)

# Order the rows by clusters, then by p-adjusted values
combined_markers <- combined_markers %>% arrange(as.numeric(as.character(cluster)), as.numeric(as.character(p_val_adj)))
combined_markers <- combined_markers %>% relocate(gene) %>% relocate(cluster)
View(combined_markers)

# Write as table
combined_out <- paste0(out_dir,"/combined_markers_",clust_no,".txt")
write.table(combined_markers,combined_out,sep="\t")

top5_comb <- combined_markers %>%
        group_by(cluster) %>%
        top_n(n = 5,
              wt = avg_log2FC)
              
# Visualize top 5 markers per cluster
View(top5_comb)

# Write as table
top5_comb_out <- paste0(out_dir,"/top5_comb_markers_",clust_no,".txt")
write.table(top5_comb,top5_comb_out,sep="\t")

# Find markers for each cluster and write to separate tables 
cell_type_list <- levels(seurat_ob)
for (cell_type in cell_type_list){
  name <- paste(cell_type,"markers",sep="_")
  name <- sub(" ", "_", name)
  name <- sub(")", "", name)
  name <- sub("-", "_", name)
  value <- FindMarkers(seurat_ob, ident.1 = cell_type)
  assign(name, value)
  out <- paste0(out_dir,"/",name,"_",clust_no,".txt")
  write.table(value,out,sep="\t")
}

### 2. Use seurat_clusters (70 clusters after doublet removal)
clust_no <- "70C"

# Set indentity classes as active clusters
Idents(object = seurat_ob) <- "active_clusters"
levels(seurat_ob)

# Find all markers
combined_markers <- FindAllMarkers(object = seurat_ob, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25)
                          
# Order the rows by clusters, then by p-adjusted values
combined_markers <- combined_markers %>%
        dplyr::arrange(as.numeric(as.character(cluster)), as.numeric(as.character(p_val_adj)))
View(combined_markers)
                          
# Write as table
combined_out <- paste0(out_dir,"/combined_markers_",clust_no)
write.table(value,combined_out,sep="\t")

# Find markers for each cluster and write to separate tables 
cell_type_list <- levels(seurat_ob)
for (cell_type in cell_type_list){
  name <- paste(cell_type,"markers",sep="_")
  name <- sub(" ", "_", name)
  name <- sub(")", "", name)
  name <- sub("-", "_", name)
  value <- FindMarkers(seurat_ob, ident.1 = cell_type)
  assign(name, value)
  out <- paste0(out_dir,"/",name,"_",clust_no)
  write.table(value,out,sep="\t")
}

# Done
