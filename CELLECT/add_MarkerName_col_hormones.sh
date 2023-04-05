#!/bin/bash

# ----------------------------------------------------------
# Add MarkerName column to hormone summary stats in preparation for mapping rsid
# MarkerName format: 1:13259:A_G
# melodyjparker14@gmail.com - Apr 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 10
#SBATCH -J .............
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

cd //well/lindgren/users/mzf347/p50/sumstats

mkdir MAF_filtered
time for f in *.txt; do awk '{if ($5 < $6) {a1=$5; a2=$6} else {a1=$6; a2=$5}; if(NR==1) {print "MarkerName",$0} else {print $2":"$3":"a1"_"a2,$0}}' $f > MAF_filtered/MN_$f; done

# When tested in interactive node this took 2m7.971s


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
