#!/bin/bash

# ----------------------------------------------------------
# Run CELLECT
# melodyjparker14@gmail.com - Apr 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH -J run_CELLECT
#SBATCH -o logs/output.out
#SBATCH -e logs/error.err
#SBATCH -a 1-6

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

IN=config

# Make an index file for the config files if it doesn't exist
if [ ! -f "$IN"/index.txt ]; then
  for f in "$IN"/*; do basename ${f} >> "$IN"/index.txt; done
fi

CONFIG_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' "$IN"/index.txt)

# Run CELLECT-LDSC
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile "$IN"/"$CONFIG_FILE"

# Run CELLECT-GENES
# snakemake --use-conda -j -s cellect-genes.snakefile --configfile "$IN"/"$CONFIG_FILE"

conda deactivate

echo "###########################################################"
echo "Slurm array task $SLURM_ARRAY_TASK_ID finished at: "`date`
echo "###########################################################"
exit 0
