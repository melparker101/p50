# CELLECT
This file goes through how to set up and use CELLEX and CELLECT and gives a few tips to avoid making the same mistakes that I did...

### 1. Ensure conda is set up
I installed miniforge/mambaforge, but BMRC recommends to load their software module.
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

- Loading BMRC module
```
# Check what versions of Anaconda are availible on the cluster
module avail Anaconda

# Load and initialise conda
# These commands would need to be run every time you want to use conda in your bash session. I would add these commands to your ~/.bashrc so that it is always automatically loaded and initiated when you start a session. 
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

- Installing mambaforge manually (BMRC does not recommend this, but it has worked fine for me)
Follow the instructions on the [mamba website](https://mamba.readthedocs.io/en/latest/installation.html) and install either mambaforge or miniforge.
```
# Add mambaforge to the END of the path (Adding it to the start of the path in .bashrc caused issues for me)
export PATH=$PATH:/users/lindgren/mzf347/mambaforge/bin/

# Restart terminal (.bashrc is automatically run when a new bash session is started)

# See if mamba has installed properly - this command should give an output
mamba

# #(make sure envs and pcks is set to not install in home dir, 
#see here under anaconda https://www.medsci.ox.ac.uk/divisional-services/support-services-1/bmrc/scientific-software-directory)
```

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
https://github.com/perslab/CELLEX#quick-start
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

### 3. Install CELLECT
IMPORTANT: Make sure to use git lfs and --recurse-submodules when cloning the CELLECT repository!
https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial
```
# Create software folder
mkdir //well/lindgren/users/yourusername/software

# Make sure git lfs is set up before you clone the repository
git lfs env 

# Clone the repository
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
```




#https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial

#go into CELLECT folder and check if ldsc folder empty, if so git lfs didn't work properly


#only run once when setting it all up
conda env create -f ldsc/environment_munge_ldsc.yml

#this you will aleays need to run when running new data
source activate activate munge_ldsc

#then run examples from tutorial -> works to munge data
#(you will have to remove dots from colnames in examples/tabula_muris-test.csv)

#to make sure second part of tutorial works run all this only once
#actually its the 3rd part, they mention CELLEX in the 2nd part but that is not set up yet but they provide example data


#make a new environment that installs snakemake

conda create -c conda-forge -c bioconda -n snakemake snakemake


pip install pybedtools 

conda install bedops

#back to tutorial

module load BEDTools/2.29.2-GCC-8.3.0

#the next command  we only run the first time on an login node for the tutorial since it
#will try to download stuff which we cant when on an interactive node
#when running with real data make sure to pack in script

snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile config.yml













2.
3.
4. p50/config_p50.yml has the following edits from the config.yml file provided by CELLECT:
BASE_OUTPUT_DIR: p50/CELLECT_OUT_p50
