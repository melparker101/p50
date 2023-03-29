# Using CELLEX on ovary datasets
Make sure conda is installed. Using conda environments is recommended for the following steps. 

1. Create this using the **seurat2h5.yml** yaml file provided using this command:
``` bash
# Create and activate conda env
conda env create --name cellex --file=seurat2h5.yml
conda activate cellex

# Load R
module load R/4.1.2-foss-2021b
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
```
Then use the corresponding R file **R_files/GSE123456.R** to prepare the count and metadata required as input for CELLEX. 
This must be done manually because the datasets are all in different formats and also need dataset-specific filtering.

2. Create this using the **cellex.yml** yaml file provided using this command:
``` bash
# Create and activate conda env
conda env create --name cellex --file=cellex.yml
conda activate cellex
```
Use **python_files/GSE123456.py** to load in the count and metadata into python then run CELLEX.
