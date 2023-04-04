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


awk 'NR==1 {print "rsid",$0; next} NR==FNR {a[$1]=$0; next} $1 in a {print $2,a[$1]}' test_ss.txt ../dbSNP/MarkerName_map_GRCh37.txt > female_infertility_UKBB.txt
# Use GRCh37_GCF_000001405.25_map.txt to get more rsids


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
