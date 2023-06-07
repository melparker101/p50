# p50 Infertility Project
### The tasks
This repository contains code used for the p50 infertility project. There are two separate tasks:
- Running [CELLECT](https://github.com/melparker101/p50_Infertility/tree/main/CELLECT) on ovary datasets with infertility and hormone GWAS sumstats to prioritise etiologic cell types.
- Finding [marker genes](https://github.com/melparker101/p50_Infertility/tree/main/cluster_marker_genes) for clusters from the ovary datasets.

### Preparing datasets
- See [datasets](https://github.com/melparker101/p50_Infertility/tree/main/datasets) for more information on the datasets used for this part of the p50 project. 
- The [create_h5seurat_h5ad](https://github.com/melparker101/p50_Infertility/tree/main/create_h5seurat_h5ad) directory contains code for creating h5seurat files (R Seurat object) from the scRNA-seq data required for finding cluster gene markers and then converting to h5ad files (python Anndata object) to prepare for CELLEX. Where cell-type annotations are not availible online, the code includes clustering and annotating cell types using Seurat. 
