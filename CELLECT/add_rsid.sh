#!/bin/bash

# ----------------------------------------------------------
# Add rsid to summary stats using an SNP map file, joining on the first column in each file: MarkerName
# melodyjparker14@gmail.com - Apr 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 10
#SBATCH -J add_rsid
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
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"


# Set working directory
# cd = //well/lindgren/users/mzf347/p50

# Set variables
IN=data/sumstats/other
OUT=data/sumstats/other
# OUT=data/sumstats/premunge
MAP=dbSNP

# Use GRCh37_GCF_000001405.25_map.txt to get more rsids

# Make an index file if it doesn't exist
if [ ! -f "$IN"/index.txt ]; then
  for f in "$IN"/MAFfiltered*; do basename ${f} >> "$IN"/index.txt; done
fi

# Make out dir if it doesn't exist
# mkdir -p "$OUT"

# Set sumstats file name for each slurm array task
SUMSTATS_FILE=$(sed "${SLURM_ARRAY_TASK_ID}"'q;d' "$IN"/index.txt)

# Add rsid column and save in new folder
awk 'NR==1 {print "RSID",$0; next} NR==FNR {a[$1]=$0; next} $1 in a {print $2,a[$1]}' \
    "$IN"/"$SUMSTATS_FILE" "$MAP"/MarkerName_map_GRCh37.txt > \
    "$OUT"/premunge_"${SUMSTATS_FILE##*MAFfiltered_MN_Ncol_}"


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
