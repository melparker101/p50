# Using CELLEX on ovary datasets
Using conda environments is recommended for the following steps; make sure conda is already installed.
(For more details about the CELLEX tool, see their [github](https://github.com/perslab/CELLEX).)

1. Create the seurat2h5 conda environment by using the **seurat2h5.yml** yaml file provided, then load R:
``` bash
# Create and activate conda env
conda env create --name seurat2h5 --file=seurat2h5.yml
conda activate seurat2h5

# Load R
module load R/4.1.2-foss-2021b
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
```
Then use the corresponding R file **R_files/GSE123456.R** to prepare the count and metadata required as input for CELLEX. 
This must be done manually because the datasets are all in different formats and also need dataset-specific filtering.
``` bash
# Deactivate conda environment when done
conda deactivate
```

2. Create the cellex conda environment using the **cellex.yml** yaml file:
``` bash
# Create and activate conda env
conda env create --name cellex --file=cellex.yml
conda activate cellex
```
Use **python_files/GSE123456.py** to load in the count and metadata into python then run CELLEX on it.
``` bash
# Deactivate conda environment when done
conda deactivate
```
