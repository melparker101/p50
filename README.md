# Using CELLEX and CELLECT on single cell RNA-seq ovary datasets with infertility GWAS summary statistics
### 1. Download datasets
We need scRNA-seq count data and the corresponding cell type annotations metadata.
### 2. Run CELLEX
Using the counts and metadata as input, we use CELLEX to produce expression specificity files (ESMU).
### 3. Run CELLECT
Using munged summary stats and ESMU files as input, we use CELLECT to prioritise etilogical cell types.
### 4. Visualisation
Use R to visualise the results.

---

For more information on each of the steps, see the README.md files of each sub-directory.
