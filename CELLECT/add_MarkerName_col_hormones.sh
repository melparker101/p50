#!/bin/bash

# ----------------------------------------------------------
# Add MarkerName column to hormone summary stats in preparation for mapping rsid
# MarkerName format: 1:13259:A_G
# melodyjparker14@gmail.com - Apr 23
# ----------------------------------------------------------

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 10
#SBATCH -J add_MN_col_hormones
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


IN="data/sumstats/other/"

# Hormomes sumstats colnames: ID CHROM GENPOS MAF Allele1 Allele2 Freq1 FreqSE BETA SE PVALUE Direction HetPVal

# Time the command
# Loop through sumstats files
# First half of awk command sorts alleles into alphabetical order and labels them a1 and a2
# Second half of awk command adds column name and then marker name concatenated strings (chr:pos:a1_a2) for the other lines

# Loop through sumstats files
time for f in "$IN"/Ncol*.txt; do
    if [ $(basename $f) = "Ncol_Infertility1_F_EUR.txt" ]; then
        cp "$IN"/Ncol_Infertility1_F_EUR.txt "$IN"/MN_Ncol_Infertility1_F_EUR.txt
        continue  # Skip adding marker name col for the infertility file - it already has one
    fi
    # Add MarkerName column to sumstats and save to new file
    awk '{if ($5 < $6) {a1=$5; a2=$6} else {a1=$6; a2=$5}; \
            if(NR==1) {print "MarkerName\t"$0} else {print $2":"$3":"a1"_"a2"\t"$0}}' $f \
                > "$IN"/MN_$(basename $f)
done

# When tested in interactive node this took 1m47.966s


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
