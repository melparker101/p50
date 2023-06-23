###########################################
# Add N column to summary stats
# melodyjparker14@gmail.com - May 23
# This script creates a "directions" file per sumstats which tells us which studies each SNP is present in.
# 1 = SNP present in study, 0 = SNP not present in study
# It then creates an N column
# It assumes there exists text files with cohort sample sizes for each hormone phenotype
# It assumes the output file exists
# This hasn't been tested.. run code chunks separately 
###########################################

###########################################
# Paths
###########################################
# Path for sumstats files
IN="data/sumstats/original"
OUT="data/sumstats/other"
COHORTS="data/sumstats/cohorts"

###########################################
# Functions
###########################################
# Function for calculating Neff from case and control sample sizes
calcNeff() {
  local n_cases=$1
  local n_controls=$2
  local N_eff=$(echo "2 / (1 / $n_cases + 1 / $n_controls)" | bc -l)
  echo "$N_eff"
}

###########################################
# Create a 'direction' file for hormone sumstats
###########################################
# Subset the hormones sumstats data to only include the 'ID' and 'Direction' columns
for f in "$IN"/*EUR_filtered.txt
  do awk '{print $1, $12}' "$f" > "$OUT/$(basename -- "${f%_filtered.txt}")_directions.txt"
done

# Loop through and edit each summary stats directions file
for phenotype in LH_F_EUR FSH_F_EUR Testosterone_F_EUR Progesterone_F_EUR Oestradiol_F_EUR Testosterone_sex_comb_EUR
do
    # Set shell variables
    COHORT1=$(awk 'NR==2 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT2=$(awk 'NR==3 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT3=$(awk 'NR==4 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT4=$(awk 'NR==5 {print $1}' "$COHORTS"/${phenotype}.txt)
    FILENAME=${phenotype}_directions.txt
    
    # Print variable names
    echo $COHORT1
    echo $COHORT2
    echo $COHORT3
    echo $COHORT4
    echo $FILENAME
    echo ""
    
    # Split up the characters of the 'Direction' column and name columns by study name
    # Use -v to pass shell variables into awk command
    # Overwrite file
    awk -v c1="$COHORT1" -v c2="$COHORT2" -v c3="$COHORT3" -v c4="$COHORT4" \
        'NR==1 {print $1, c1, c2, c3, c4; next} NR==FNR {printf "%s ", $1; gsub(/.{1}/,"& ",$2); print $2}' "$OUT"/"$FILENAME" \
            > tmp && mv tmp "$OUT"/"$FILENAME"

    # Change ? to 0 and +/- to 1
    awk '{gsub(/\?/,"0"); gsub(/[+-]/,"1"); print}' "$OUT"/"$FILENAME" > tmp && mv tmp "$OUT"/"$FILENAME"
done

###########################################
# Add N column to hormone sumstats
###########################################
# Create a sym link to link the infertility file to the same file name format as the hormones sumstats
ln -s $(readlink -fe "$IN"/female_infertility_analysis1_eur_MA_results_chr_pos.txt) "$IN"/Infertility1_F_EUR_filtered.txt
ln -s $(readlink -fe "$IN"/female_infertility_analysis2_eur_MA_results_chr_pos.txt) "$IN"/Infertility2_F_EUR_filtered.txt
ln -s $(readlink -fe "$IN"/female_infertility_analysis3_eur_MA_results_chr_pos.txt) "$IN"/Infertility3_F_EUR_filtered.txt
ln -s $(readlink -fe "$IN"/female_infertility_analysis4_eur_MA_results_chr_pos.txt) "$IN"/Infertility4_F_EUR_filtered.txt
ln -s $(readlink -fe "$IN"/female_infertility_analysis5_eur_MA_results_chr_pos.txt) "$IN"/Infertility5_F_EUR_filtered.txt

# Add N column to sumstats and save to new file - Just do for hormones for now
for phenotype in LH_F_EUR FSH_F_EUR Testosterone_F_EUR Progesterone_F_EUR Oestradiol_F_EUR Testosterone_sex_comb_EUR # Infertility1_F_EUR
do
    COHORT1=$(awk 'NR==2 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT2=$(awk 'NR==3 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT3=$(awk 'NR==4 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT4=$(awk 'NR==5 {print $1}' "$COHORTS"/${phenotype}.txt)
    SAMPLE_SIZE1=$(awk 'NR==2 {print $2}' "$COHORTS"/${phenotype}.txt)
    SAMPLE_SIZE2=$(awk 'NR==3 {print $2}' "$COHORTS"/${phenotype}.txt)
    SAMPLE_SIZE3=$(awk 'NR==4 {print $2}' "$COHORTS"/${phenotype}.txt)
    SAMPLE_SIZE4=$(awk 'NR==5 {print $2}' "$COHORTS"/${phenotype}.txt)
    FILENAME=${phenotype}_directions.txt
    
    echo ${phenotype}
    echo "Sample size 1:${SAMPLE_SIZE1}" 
    echo "Sample size 2:${SAMPLE_SIZE2}"
    echo "Sample size 3:${SAMPLE_SIZE3}"
    echo "Sample size 4:${SAMPLE_SIZE4}"
  
    paste "$IN"/${phenotype}_filtered.txt \
        <(awk -v n1="$SAMPLE_SIZE1" -v n2="$SAMPLE_SIZE2" -v n3="$SAMPLE_SIZE3" -v n4="$SAMPLE_SIZE4" \
            'NR==1 {print "N"; next} NR>1 { print $2*n1 + $3*n2 + $4*n3 + $5*n4 }' "$OUT"/${phenotype}_directions.txt) \
                > "$OUT"/Ncol_${phenotype}.txt
done

###########################################
# Add N column to infertility sumstats
###########################################
# Function for calculating Neff from case and control sample sizes
calcNeff() {
  local n_cases=$1
  local n_controls=$2
  local N_eff=$(echo "2 / (1 / $n_cases + 1 / $n_controls)" | bc -l)
  echo "$N_eff"
}

for phenotype in Infertility1_F_EUR Infertility2_F_EUR Infertility3_F_EUR Infertility4_F_EUR Infertility5_F_EUR
do
  #awk '{print $0, calcNeff $9 $10}' "$IN"/${phenotype}_filtered.txt | head -3
  #awk '{print $0, calcNeff($9, $10)}' "$IN/${phenotype}_filtered.txt" | head -3
  paste "$IN"/${phenotype}_filtered.txt \
      <(awk '{print $0, system(calcNeff $9 $10)}' "$IN/${phenotype}_filtered.txt") | head -3
done

for phenotype in Infertility1_F_EUR Infertility2_F_EUR Infertility3_F_EUR Infertility4_F_EUR Infertility5_F_EUR
do
  paste "$IN/${phenotype}_filtered.txt" \
      <(awk -v phenotype="${phenotype}" '{print $0, calcNeff($9, $10)}' "$IN/${phenotype}_filtered.txt" | head -3) \
      > "$OUT/Ncol_${phenotype}.txt"
done
##########################################

# Function for calculating Neff from case and control sample sizes
export -f calcNeff

for phenotype in Infertility1_F_EUR Infertility2_F_EUR Infertility3_F_EUR Infertility4_F_EUR Infertility5_F_EUR
do
  paste "$IN/${phenotype}_filtered.txt" \
      <(awk -v phenotype="${phenotype}" '{print $0, "bash -c '\''calcNeff $9 $10'\''"}' "$IN/${phenotype}_filtered.txt" | bash ) | head -3
done


  

###########################################
# End
###########################################
