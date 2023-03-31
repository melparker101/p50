# Creating Conda Environments
After creating the conda environments, check that all of the important packages have properly installed. 

### 1. Download Anaconda or Miniconda
- For the BMRC cluster, see the [staff resources](https://www.medsci.ox.ac.uk/for-staff/resources/bmrc/python-on-the-bmrc-cluster).
- Otherwise, see the [conda user guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

### 2. Create the conda environments
- CELLEX
Use the **cellex.yml** yaml file provided here to create the cellex conda environment. 
This environment includes packages required for running CELLEX and also includes scanpy which is used for data formatting in python.
``` bash
# Create and activate conda env
conda env create --name cellex --file=cellex.yml
conda activate cellex
```

- Sumstats Munging
Create a munge_ldsc environment which uses python2 using [environment.yml](https://github.com/pascaltimshel/ldsc/blob/d869cfd1e9fe1abc03b65c00b8a672bd530d0617/environment.yml)
``` bash
conda env create -f ldsc/environment_munge_ldsc.yml # creates 'munge_ldsc' environment 
```

- CELLECT
A snakemake environment is required to run CELLECT. See [CELLECT](https://github.com/perslab/CELLECT) for more information.
``` bash
conda install -c bioconda -c conda-forge snakemake">=5.27.4"
```
