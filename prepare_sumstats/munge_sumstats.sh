#!/bin/bash

# ----------------------------------------------------------
# Script to munge sumstats
# Run inside p50
# //well/lindgren/users/mzf347/CELLECT/p50
# The CELLECT directory is cloned from their github
# melodyjparker14@gmail.com - Apr 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J munge_sumstats
#SBATCH -o logs/output.out
#SBATCH -e logs/error.err
#SBATCH -a 1-7

#  Parallel environment settings 
#  For more information on these please see the documentation 
#  Allowed parameters: 
#   -c, --cpus-per-task 
#   -N, --nodes 
#   -n, --ntasks 

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

# Activate munge conda env.
source ~/.bashrc
conda activate munge_ldsc

# Define paths
CELLECT='..'  # //well/lindgren/users/mzf347/CELLECT
IN=data/sumstats/other
OUT=data/sumstats/munged

# Make an index file for munging
if [ ! -f "$IN"/munge_index.txt ]; then
  for f in "$IN"/premunge*; do basename ${f##*premunge_} >> "$IN"/munge_index.txt; done
fi

SUMSTATS_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' "$IN"/munge_index.txt)

# Make outdir
mkdir -p "$OUT"

# Run munge script
python "$CELLECT"/ldsc/mtag_munge.py \
--sumstats "$IN"/premunge_"$SUMSTATS_FILE" \
--merge-alleles "$CELLECT"/data/ldsc/w_hm3.snplist \
--keep-pval \
--ignore MarkerName \
--out "$OUT"/munged_"${SUMSTATS_FILE%.txt}"

# Deactivate conda environment
conda deactivate

echo "###########################################################"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
