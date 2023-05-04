# Set Up
Download the required packages and create the recommended conda environments.

## 1. Load or download Anaconda or Miniconda
- For the BMRC cluster, see the [staff resources](https://www.medsci.ox.ac.uk/for-staff/resources/bmrc/python-on-the-bmrc-cluster).
- See the [conda user guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

## 2. Create conda environments
### CELLEX

Install [CELLEX](https://github.com/perslab/CELLEX) and [scanpy](https://scanpy.readthedocs.io/en/stable/installation.html).
```bash
pip install cellex
pip install scanpy
```
If this doesn't work, create the following cellex environment using [cellex.yml](https://github.com/melparker101/p50/blob/main/set_up/cellex.yml).
CELLEX ran successfully using in this environment with the python scripts in this repository.
``` bash
# Create and activate conda env
conda env create --name cellex --file=cellex.yml
conda activate cellex
```

### Sumstats munging
Create a munge_ldsc environment which uses python2 using [environment.yml](https://github.com/pascaltimshel/ldsc/blob/d869cfd1e9fe1abc03b65c00b8a672bd530d0617/environment.yml).
``` bash
conda env create -f ldsc/environment_munge_ldsc.yml  # creates 'munge_ldsc' environment 
```

### CELLECT
A snakemake environment is required to run [CELLECT](https://github.com/perslab/CELLECT).
``` bash
conda install -c bioconda -c conda-forge snakemake">=5.27.4"
```
