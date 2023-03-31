# Using CELLEX and CELLECT on single cell RNA-seq ovary datasets with infertility GWAS summary statistics
### 1. Download datasets
We need scRNA-seq count data and the corresponding cell type annotations metadata. See [p50/datasets](https://github.com/melparker101/p50/tree/main/datasets) for more information.
### 2. Run CELLEX
Using the counts and metadata as input, we use [CELLEX](https://github.com/perslab/CELLEX) to produce expression specificity files (ESMU). See [p50/CELLEX](https://github.com/melparker101/p50/tree/main/CELLEX) for R and python code used to prepare data and run CELLEX.
### 3. Run CELLECT
Using munged summary stats and ESMU files as input, we use [CELLECT](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial) to prioritise etilogical cell types. See [p50/CELLECT](https://github.com/melparker101/p50/tree/main/CELLECT).
### 4. Visualisation
Use R to visualise the results. See [p50/visualisation](https://github.com/melparker101/p50/tree/main/visualisation).

---

For more information on each of these steps, see the README.md files of each sub-directory.
