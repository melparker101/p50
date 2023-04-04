#!/bin/bash

# ----------------------------------------------------------
# Manipulate map file ready for use on summary stats
# This has been tested and works
# melodyjparker14@gmail.com - Apr 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 10
#SBATCH -J manipulate_map_file
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


# 1. filter out MT, Y, N
sed -e '/H/d' -e '/MT/d' -e '/Y/d' GRCh37_GCF_000001405.25_map.txt > 1-23_filtered_map.txt


# 2. filter out multi-alleles: filter out any fields that are more than one using , and " " as field separators
awk -F" |," 'length($4) <= 1 && length($5) <= 1 && length($6) <= 1 && length($7) <= 1' 1-23_filtered_map.txt > sa_1-23_filtered_map.txt


# 3. Print out a new line for each SNP combination for each rsid
awk -F" |," '{
if(!$6)
{	print $0
}
else if(!$7)
{	print $1,$2,$3,$4,$5;
	print $1,$2,$3,$4,$6
}
else 
{	print $1,$2,$3,$4,$5;
	print $1,$2,$3,$4,$6;
	print $1,$2,$3,$4,$7
}
}' sa_1-23_filtered_map.txt > final_map.txt

# 4. Make into a marker name file
awk '{if ($4 < $5) {a1=$4; a2=$5} else {a1=$5; a2=$4}; print $1":"$2":"a1"_"a2,$3}' final_map.txt > MarkerName_map_GRCh37.txt


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
