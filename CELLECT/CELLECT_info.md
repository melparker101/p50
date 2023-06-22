# CELLECT Tutorial

## CELLEX and CELLECT Set Up
This markdown file goes through how to set up and use CELLEX and CELLECT and gives a few tips to avoid making the same mistakes that I did...

### General tips:
- **Run CELLECT on rescomp the first time** as it creates conda environments (internet connection is required for necessary package installation). After this, send slurm scripts off. CELLECT creates these conda environments in a directory called '.snakemake'.
- For CELLECT-ldsc, make sure the ldsc folder is properly downloaded - 'git clone' does not download this directory. **Make sure --recurse-submodules argument is used** to clone the ldsc directory.
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
- For CELLECT-magma, use the 'keep p-val' option for munging - the summary stats must have a p-val column
- When running CELLECT-magma, I was given an error message "ModuleNotFoundError: No module named 'statsmodels'; import statsmodels.api as sm". Although CELLECT is meant to install all packages automatically on the first run (apart from snakemake), I had to install this manually (conda install statsmodels) in the snakemake conda environment before attempting to run CELLECT-magma again - this fixed the issue.
- For CELLECT-genes, use this temporary bug fix: perslab/CELLECT#81
- Make sure summary stats do not contain "." in the column names (there was an error caused by this in the tutorial when I tried it)
- When preprocessing the count data in R and saving as a h5ad file ready for python ready for input for CELLEX, remove the scale data from the Seurat object first. There is an issue with Seurat disk and the raw counts do not get transferred to the h5ad if not! See https://github.com/mojaveazure/seurat-disk/issues/75 for more info.

### 1. Ensure conda is set up
I installed miniforge/mambaforge manually, but BMRC recommends not to do this and to load their software module instead (it has worked for me though).
https://www.medsci.ox.ac.uk/for-staff/resources/bmrc/python-on-the-bmrc-cluster

By default, conda stores environments and packages in your home directory, but we do not have much space here. We want to change this configuration and redirect this.
```
# Create a directory to store conda enviromnents
mkdir /well/group/users/username/conda

# Create a .condarc file (conda configuration file)
cat > ~/.condarc

channels:
  - conda-forge
  - bioconda
  - defaults
 
pkgs_dirs:
  - /well/group/users/username/conda/${MODULE_CPU_TYPE}/pkgs
envs_dirs:
  - /well/group/users/username/conda/${MODULE_CPU_TYPE}/envs
  
# ctrl D  
```

**(i) Loading BMRC module**
```
# Check what versions of Anaconda are availible on the cluster
module avail Anaconda

# Load and initialise conda
# These commands would need to be run every time you want to use conda in your bash session.
# I would add these commands to your ~/.bashrc so that it is always automatically loaded and initiated when you start a session. 
module load Anaconda3/2022.05
eval "$(conda shell.bash hook)"
```

When using BMRC, make sure you make identical conda environments for both skylake (cluster1 or cluster2) and ivybridge (rescomp3).
```
# To access rescomp 3
ssh rescomp3
```

```
module use -a /apps/eb/dev/{skylake,ivybridge}/modules/all
# By executing this command, you are adding the specified module paths (/apps/eb/dev/skylake/modules/all and /apps/eb/dev/ivybridge/modules/all) to the module search path.
```

**(ii) Installing mambaforge manually (BMRC does not recommend this, but it worked for me)**
Follow the instructions on the [mamba website](https://mamba.readthedocs.io/en/latest/installation.html) and install either mambaforge or miniforge.
```
# Add mambaforge to the END of the path (Adding it to the start of the path in .bashrc caused issues for me)
export PATH=$PATH:/users/lindgren/mzf347/mambaforge/bin/

# Restart terminal (.bashrc is automatically run when a new bash session is started)

# See if mamba has installed properly - this command should give an output
mamba
```
Make sure envs and pcks is set to not install in home dir, see above and [here](https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/scientific-software-directory).

My ~/.bashrc contains the following lines which initiate conda automatically when I start my bash session.
```
export PATH=$PATH:/users/lindgren/mzf347/mambaforge/bin/

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/gpfs3/users/lindgren/mzf347/mambaforge/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/gpfs3/users/lindgren/mzf347/mambaforge/etc/profile.d/conda.sh" ]; then
        . "/gpfs3/users/lindgren/mzf347/mambaforge/etc/profile.d/conda.sh"
    else
        export PATH="/gpfs3/users/lindgren/mzf347/mambaforge/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/gpfs3/users/lindgren/mzf347/mambaforge/etc/profile.d/mamba.sh" ]; then
    . "/gpfs3/users/lindgren/mzf347/mambaforge/etc/profile.d/mamba.sh"
fi
```

### 2. Install CELLEX
Read the [CELLEX github](https://github.com/perslab/CELLEX#quick-start).
I created a conda environment and installed CELLEX inside it.
```
# Create conda env for cellex
conda create -n cellex
conda activate cellex

# Install pip (https://askubuntu.com/questions/1025189/pip-is-not-working-importerror-no-module-named-pip-internal)
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3 get-pip.py --user --force-reinstall
pip -V

# Install CELLEX
pip install cellex
```
I can't specifically remember the issues I came across when doing this, but I know there are certain [requirements](https://github.com/perslab/CELLEX/blob/master/requirements.txt) for python package versions:
```
setuptools>=41.2.0
setuptools_scm>=3.3.3
numpy<1.24
scipy>=1.3.1
pandas>=1.0.3
h5py>=2.9.0
loompy>=3.0.6
adjustText>=0.7.3
plotnine>=0.6.0
```
Follow the CELLEX tutorials provided online: [demo_mousebrain_vascular_cells](https://github.com/perslab/CELLEX/blob/master/tutorials/demo_mousebrain_vascular_cells.ipynb) and [demo_moca_100k](https://github.com/perslab/CELLEX/blob/master/tutorials/demo_moca_100k.ipynb). We will need the ESMU output files that CELLEX produces to run the CELLECT example.

### 3. Install CELLECT
**IMPORTANT: Make sure to use git lfs and --recurse-submodules when cloning the CELLECT repository!**
Have a read of the [CELLECT-LDSC tutorial](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial). Before running on our own data, set up and run their example on rescomp. This will download some extra conda environments required to run CELLECT. Once they are downloaded, CELLECT can be used offline (we can run it using slurm scripts).

Clone the CELLECT github.
```
# Create software folder
mkdir //well/lindgren/users/yourusername/software

# Make sure git lfs is set up before you clone the repository
git lfs env 

# Clone the repository
git clone --recurse-submodules https://github.com/perslab/CELLECT.git

# Check that the ldsc folder isn't empty, if so then git lfs didn't work properly
```

Create a munging conda environment. This environment is for running their munging script on the summary statistics in, which requires python 2.
```
# Create munge env
conda env create -f ldsc/environment_munge_ldsc.yml
# The environment will be called "munge_ldsc"
```
Now make another conda enviroment containing snakemake for running CELLECT in. We can then run CELLECT using snakemake. **IMPORTANT: AN INTERNET CONNECTION IS REQUIRED THE FIRST TIME YOU RUN CELLECT - run their example on rescomp when running for the first time**
```
# Make a new environment with snakemake
conda create -c conda-forge -c bioconda -n snakemake snakemake

conda activate snakemake

# Install other packages?
pip install pybedtools 
conda install bedops
module load BEDTools/2.29.2-GCC-8.3.0

# Run their example on rescomp 
# IMPORTANT: AN INTERNET CONNECTION IS NEEDED THE FIRST TIME YOU RUN IT
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile config.yml
```
## Running CELLEX (preparing ESMU files)
Now that the CELLEX and CELLECT are set up and the example given by the authors has run smoothly, it is time to use them for our own choice of datasets. The next step is running CELLEX to generate ESMU files.
CELLEX requires as input:
- Raw count data
- Cell type annotation data

### 1. Filter and prepare count and cell-type annotations data
I used Seurat in R to filter and manipulate the count data beforehand. If cell-type annotations were not provided as metadata, I performed the clustering and annotations in Seurat. I saved the annotated Seurat object as a H5Seurat file, then converted this to a h5ad file where the data can be easily extracted using scanpy in python. I found that using the R code below was best for reproducability when saving the Seurat object and converting to h5ad format:
```
# Remove scale data (there is an issue with seurat disk and the raw counts when converting to h5ad if not) - https://github.com/mojaveazure/seurat-disk/issues/75
seurat_ob <- DietSeurat(seurat_ob)

# Make sure any metadata that is of type factor is converted to charater
# https://github.com/mojaveazure/seurat-disk/issues/23
i <- sapply(seurat_ob@meta.data, is.factor)
seurat_ob@meta.data[i] <- lapply(seurat_ob@meta.data[i], as.character)

# Save as a h5ad file
SaveH5Seurat(seurat_ob, filename = out_seurat)
Convert(out_seurat, dest = "h5ad")
```
### 2. Run CELLEX
1. Activate the cellex conda env - also make sure scanpy is installed.
2. Load the h5ad file into python with scanpy.
3. Extract the count data and cell type annotations metadata from the scanpy AnnData object (the python version of a Seurat object in R).
4. Run CELLEX.
5. Save the ESMU file.

The python scripts I used are deposited [here](https://github.com/melparker101/p50_Infertility/tree/main/CELLECT/prepare_esmu).

## Preparing the summary statistics
The summary stats need to be in a specific format to use as input for CELLECT. The
[CELLECT Input-&-Output](https://github.com/perslab/CELLECT/wiki/Input-&-Output) page states which columns are required:
```
For CELLECT-LDSC, the required columns are:

SNP: the unique SNP identifier (e.g. rsID number)
N: sample size (which may vary from SNP to SNP)
Z: the Z-score associated with the SNP effect sizes for the GWAS Additional columns are allowed but will be ignored.
For CELLECT-MAGMA, the required columns are:

SNP: the unique SNP identifier (e.g. rsID number)
N: sample size (which may vary from SNP to SNP)
PVAL: the P-value associated with the SNP effect sizes for the GWAS Additional columns are allowed but will be ignored.
```
**IMPORTANT: the SNP column must contain rsids**. Follow this workflow to prepare the sumstats:

1. Add N column to sumstats ([add_N_col.sh](https://github.com/melparker101/p50/blob/main/CELLECT/prepare_sumstats/add_N_col.sh))
2. Create a map file for chr:pos:a1_a2 to rsid (dbSNP/MarkerName_map_GRCh37.txt) ([create_snp_map.md](https://github.com/melparker101/p50/blob/main/CELLECT/prepare_sumstats/create_snp_map.md) and [manipulate_map_file.sh](https://github.com/melparker101/p50/blob/main/CELLECT/prepare_sumstats/manipulate_map_file.sh))
2. Add MarkerName column to hormones sumstats files (infertility sumstats already contains this) ([add_MarkerName_col_hormones.sh](https://github.com/melparker101/p50/blob/main/CELLECT/prepare_sumstats/add_MarkerName_col_hormones.sh))
3. Filter for MAF >1% in R ([MAF_filter_all_sumstats.R](https://github.com/melparker101/p50/blob/main/CELLECT/prepare_sumstats/MAF_filter_all_sumstats.R))
4. Add rsid column to all sumstats using map file ([add_rsid.sh](https://github.com/melparker101/p50/blob/main/CELLECT/prepare_sumstats/add_rsid.sh))
5. Munge using ldsc munge script ([munge_sumstats.sh](https://github.com/melparker101/p50/blob/main/CELLECT/prepare_sumstats/munge_sumstats.sh))

When munging the sumstats, make sure you use the argument to keep PVAL - even though this sumstats column is not required for CELLECT-ldsc, it is for CELLECT-magma and thus for CELLECT-genes. 

## Preparing summary statistics files
## Run CELLECT
1. I have edited the config.yml file provided by CELLECT to include the ovary datasets and summary statistics [p50/config_p50.yml](https://github.com/melparker101/p50_Infertility/blob/main/CELLECT/run_cellect/config_p50.yml) in three places:
```
BASE_OUTPUT_DIR: p50/CELLECT_OUT_p50
```
```
SPECIFICITY_INPUT:
  - id: GSE118127
    path: p50/data/esmu/GSE118127_ens.esmu.csv.gz
  - id: GSE202601
    path: p50/data/esmu/GSE202601_ens.esmu.csv.gz
  - id: GSE213216
    path: p50/data/esmu/GSE213216_ens.esmu.csv.gz
  - id: GSE189960
    path: p50/data/esmu/GSE189960_ens.esmu.csv.gz
  - id: GSE192722
    path: p50/data/esmu/GSE192722_ens.esmu.csv.gz
  - id: GSE206143
    path: p50/data/esmu/GSE206143_ens.esmu.csv.gz
```
```
GWAS_SUMSTATS:
  - id: Testosterone_F_EUR
    path: p50/data/sumstats/munged/munged_Testosterone_F_EUR.sumstats.gz
  - id: Testosterone_sex_comb_EUR
    path: p50/data/sumstats/munged/munged_Testosterone_sex_comb_EUR.sumstats.gz
  - id: LH_F_EUR
    path: p50/data/sumstats/munged/munged_LH_F_EUR.sumstats.gz
  - id: FSH_F_EUR
    path: p50/data/sumstats/munged/munged_FSH_F_EUR.sumstats.gz
  - id: Infertility1_F_EUR
    path: p50/data/sumstats/munged/munged_Infertility1_F_EUR.sumstats.gz
  - id: Progesterone_F_EUR
    path: p50/data/sumstats/munged/munged_Progesterone_F_EUR.sumstats.gz
  - id: Oestradiol_F_EUR
    path: p50/data/sumstats/munged/munged_Oestradiol_F_EUR.sumstats.gz
```
2. Go back into the CELLECT directory (/p50/..) and run the [run_CELLECT.sh](https://github.com/melparker101/p50_Infertility/blob/main/CELLECT/run_cellect/run_CELLECT.sh) script on the cluster so that we can use more nodes. This runs
- CELLECT-ldsc
- CELLECT-magma
- CELLECT-genes
using the same config file. 
MAKE SURE THE EXAMPLES OF ALL VERSIONS OF CELLECT HAVE BEEN SUCESSFULLY RUN ON RESCOMP WHERE THERE IS INTERNET FIRST - the first time it is run it creates conda environments that are required. There is less data in the example files so it will take less time to run than all of our data which we will send off to slurm.
3. Visualise the results in R or python. Plot the log of the p-values.

