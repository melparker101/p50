# Prepare ESMU files
Use [CELLEX](https://github.com/perslab/CELLEX) to generate the expression specificity (ESMU) files required as input for CELLECT. This must be done manually because the datasets are all in different formats and also need dataset-specific filtering. One of the datasets did not need this step.

### 1. Use the R code to format and filter the data, if necessary
First, preprocess the count data using the R scripts, if needed - save in H5AD format. The following versions of R and Bioconductor were used:
``` bash
# Load R
module load R/4.1.2-foss-2021b
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
Rscript <file_name.R>
```

### 2. Use python to run CELLEX
Second, use the python scripts to run CELLEX on the datasets. Make sure the correct python packages are installed (scanpy and CELLEX) or use the conda environment (see [set_up](https://github.com/melparker101/p50/blob/main/set_up)). Use the corresponding python scripts to load in the h5ad files, extract the count and metadata and then run CELLEX on it.

``` bash
conda activate cellex
sh <file_name.py>
conda deactivate
```
