#!/bin/bash

# ----------------------------------------------------------
# Script to munge sumstats
# Run inside CELLECT DIRECTORY
# //well/lindgren/users/mzf347/p50/CELLECT
# Try from p50 and state cellex path
# melodyjparker14@gmail.com - Apr 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J munge_sumstats
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
IN=//well/lindgren/users/mzf347/p50/sumstats/premunge
OUT=//well/lindgren/users/mzf347/p50/sumstats/munged
CELLECT=//well/lindgren/users/mzf347/p50/CELLECT

# Use sample_size.txt as an index file. col1 = text files, col2 = max sample sizes
# -v is to pass an external shell variables to an awk; NR is for line number; {print $j} is to print column j
SUMSTATS_FILE=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i {print $1}' $IN/sample_size.txt)
SAMPLE_SIZE=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i {print $2}' "$IN"/sample_size.txt)

# Make outdir
mkdir -p "$OUT"

# Run munge script
python "$CELLECT"/ldsc/mtag_munge.py \
--sumstats "$IN"/"$SUMSTATS_FILE" \
--merge-alleles "$CELLECT"/data/ldsc/w_hm3.snplist \
--n-value "$SAMPLE_SIZE" \
--ignore MarkerName \
--out "$OUT"/munged_FSH_F_EUR

echo "###########################################################"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
