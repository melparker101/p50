# Using CELLEX on ovary datasets
The CELLEX conda environment is recommended for this step. Make sure conda is installed, then create this using the yaml file provided using this command:
``` bash
conda env create --name cellex --file=cellex.yml
conda activate cellex
```
1. Use the corresponding R file R_files/accession.R to prepare the count and metadata required as input for CELLEX. 
This must be done manually because the datasets are all in different formats and also need dataset-specific filtering.
2. Use python_files/accession.py to load in the count and metadata into python then run CELLEX.
