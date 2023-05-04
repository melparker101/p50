# p50
This repository contains code used for the p50 infertility project. There are two separate tasks:
- Running CELLECT on ovary datasets with infertility and hormone GWAS sumstats to prioritise etiologic cell types.
- Finding marker genes for clusters from the ovary datasets.

## Using CELLEX and CELLECT on single cell RNA-seq ovary datasets with infertility GWAS summary statistics
### 1. Download datasets
We need scRNA-seq count data and the corresponding cell type annotations metadata. See [datasets](https://github.com/melparker101/p50/tree/main/datasets) for more information.
### 2. Set up environments
Download required packages and create the recommended conda environments. More information is given in [set_up](https://github.com/melparker101/p50/tree/main/set_up).
### 3. Run CELLEX
Using the counts and metadata as input, we use [CELLEX](https://github.com/perslab/CELLEX) to produce expression specificity files (ESMU). See [CELLEX](https://github.com/melparker101/p50/tree/main/CELLEX) for R and python code used to prepare data and run CELLEX.
### 4. Run CELLECT
Using munged summary stats and ESMU files as input, we use [CELLECT](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial) to prioritise etilogical cell types. See [CELLECT](https://github.com/melparker101/p50/tree/main/CELLECT).
### 5. Visualisation
Use R to visualise the results. See [visualisation](https://github.com/melparker101/p50/tree/main/visualisation).

---

## Finding marker genes for clusters
See [marker_genes](https://github.com/melparker101/p50/tree/main/marker_genes) scripts.

## Directory structure
The basic directory stucture:
```
tree p50
|-- CELLECT_OUT_p50
|   |-- CELLECT-GENES
|   |-- CELLECT-LDSC
|   `-- CELLECT-MAGMA
|-- cluster_markers
|   |-- GSE118127
|   |-- GSE202601
|   `-- GSE213216
|-- data
|   |-- counts
|   |-- esmu
|   `-- sumstats
|-- dbSNP
|-- logs
`-- plots
```

