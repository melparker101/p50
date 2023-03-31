# Running CELLEX on ovary datasets

# 1. Use the R scripts to format and filter the data, if necessary
The following versions of R and Bioconductor were used:
``` bash
# Load R
module load R/4.1.2-foss-2021b
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
```

2. Use python to run CELLEX
Make sure the correct packages are installed or use the conda environment (see [set_up](https://github.com/melparker101/p50/blob/main/set_up).

Then use the corresponding R file to prepare the count and metadata required as input for CELLEX. 
This must be done manually because the datasets are all in different formats and also need dataset-specific filtering. One of the datasets did not need this step.

Use the corresponding python scripts to load in the h5ad files, extract the count and metadata and then run CELLEX on it.

``` bash
conda activate cellex

# Run code from appropriate python file

conda deactivate
```
