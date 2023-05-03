#!/bin/bash

# ----------------------------------------------------------
# Run CELLECT
# melodyjparker14@gmail.com - Apr 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 12
#SBATCH -J run_CELLECT
#SBATCH -o logs/CELLECT.out
#SBATCH -e logs/CELLECT.err

#  Parallel environment settings
#  For more information on these please see the documentation
#  Allowed parameters:
#   -c, --cpus-per-task
#   -N, --nodes
#   -n, --ntasks

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "##########################################################"

source ~/.bashrc
conda activate snakemake

echo $CONDA_DEFAULT_ENV

IN=p50
CONFIG_FILE=config_p50.yml

# Run CELLECT-LDSC
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile "$IN"/"$CONFIG_FILE"

# Run CELLECT-MAGMA
snakemake --use-conda -j -s cellect-magma.snakefile --configfile "$IN"/"$CONFIG_FILE"

# Run CELLECT-GENES
snakemake --use-conda -j -s cellect-genes.snakefile --configfile "$IN"/"$CONFIG_FILE"

conda deactivate

echo "###########################################################"
echo "Slurm array task $SLURM_ARRAY_TASK_ID finished at: "`date`
echo "###########################################################"
exit 0
