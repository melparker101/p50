###########################################
# Add N column to summary stats
# melodyjparker14@gmail.com - May 23
# This script creates a "directions" file per sumstats which tells us which studies each SNP is present in.
# 1 = SNP present in study, 0 = SNP not present in study
# It then creates an N column
# It assumes there exists text files with cohort sample sizes for each hormone phenotype
###########################################

###########################################
# Paths
###########################################
# Path for sumstats files
IN=data/sumstats/original
OUT=data/sumstats/other
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
# Create a 'direction' file
###########################################
# Subset the hormones sumstats data to only include the 'ID' and 'Direction' columns
for f in "$IN"/*EUR_filtered.txt
  do awk '{print $1, $12}' "$f" > "$OUT/$(basename -- "${f%_filtered.txt}")_directions.txt"
done

# Subset the infertility sumstats data to only include the 'MarkerName' and 'Direction' columns
awk '{print $1, $11}' "$IN"/female_infertility_analysis1_UKBB_Finngen_EstBB_noMACfilter_March20231.out > "$OUT"/Infertility1_F_EUR_directions.txt

# Loop through and edit each summary stats directions file
for phenotype in LH_F_EUR FSH_F_EUR Testosterone_F_EUR Progesterone_F_EUR Oestradiol_F_EUR Testosterone_sex_comb_EUR Infertility1_F_EUR
do
    # Set shell variables
    COHORT1=$(awk 'NR==2 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT2=$(awk 'NR==3 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT3=$(awk 'NR==4 {print $1}' "$COHORTS"/${phenotype}.txt)
    FILENAME=${phenotype}_directions.txt
    
    # Print variable names
    echo "$COHORT1 $COHORT2 $COHORT3 $FILENAME"
    
    # Split up the characters of the 'Direction' column and name columns by study name
    # Use -v to pass shell variables into awk command
    # Overwrite file
    awk -v c1="$COHORT1" -v c2="$COHORT2" -v c3="$COHORT3" \
        'NR==1 {print $1, c1, c2, c3; next} NR==FNR {printf "%s ", $1; gsub(/.{1}/,"& ",$2); print $2}' "$OUT"/"$FILENAME" \
            > tmp && mv tmp "$OUT"/"$FILENAME"

    # Change ? to 0 and +/- to 1
    awk '{gsub(/\?/,"0"); gsub(/[+-]/,"1"); print}' "$OUT"/"$FILENAME" > tmp && mv tmp "$OUT"/"$FILENAME"
done

###########################################
# Calculate Neff for infertility sumstats
###########################################
# Infertility1 case/control variables:
# Define FinnGen variables manually
FinnGen_N_Cases=14759
FinnGen_N_Controls=111583
FinnGen_N_eff=$(calcNeff $FinnGen_N_Cases $FinnGen_N_Controls)  # 26069.77

# Define UKBB variables
UKBB_N_Cases=3447
UKBB_N_Controls=245591
UKBB_N_eff=$(calcNeff $UKBB_N_Cases $UKBB_N_Controls)  # 6798.578

# Define EstBB variables
EstBB_N_Cases=12027
EstBB_N_Controls=111395
EstBB_N_eff=$(calcNeff $EstBB_N_Cases $EstBB_N_Controls)  # 21710.03

# Create an Neff sample size file for infertility
echo "cohort sample_size" > $COHORTS/Infertility1_F_EUR.txt
echo "FinnGen" $FinnGen_N_eff >> $COHORTS/Infertility1_F_EUR.txt
echo "UKBB" $UKBB_N_eff >> $COHORTS/Infertility1_F_EUR.txt
echo "EstBB" $EstBB_N_eff >> $COHORTS/Infertility1_F_EUR.txt

###########################################
# Add N column to sumstats
###########################################
# Create a sym link to link the infertility file to the same file name format as the hormones sumstats
ln -s $(readlink -fe "$IN"/female_infertility_analysis1_UKBB_Finngen_EstBB_noMACfilter_March20231.out) "$IN"/Infertility1_F_EUR_filtered.txt

# Add N column to sumstats and save to new file
for phenotype in LH_F_EUR FSH_F_EUR Testosterone_F_EUR Progesterone_F_EUR Oestradiol_F_EUR Testosterone_sex_comb_EUR Infertility1_F_EUR
do
    COHORT1=$(awk 'NR==2 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT2=$(awk 'NR==3 {print $1}' "$COHORTS"/${phenotype}.txt)
    COHORT3=$(awk 'NR==4 {print $1}' "$COHORTS"/${phenotype}.txt)
    SAMPLE_SIZE1=$(awk 'NR==2 {print $2}' "$COHORTS"/${phenotype}.txt)
    SAMPLE_SIZE2=$(awk 'NR==3 {print $2}' "$COHORTS"/${phenotype}.txt)
    SAMPLE_SIZE3=$(awk 'NR==4 {print $2}' "$COHORTS"/${phenotype}.txt)
    FILENAME=${phenotype}_directions.txt
    
    # awk -v n1="$SAMPLE_SIZE1" -v n2="$SAMPLE_SIZE2" -v n3="$SAMPLE_SIZE3" '{ print $1*n1 + $2*n2 + $3*n3 }' "$OUT"/${phenotype}_directions.txt > "$OUT"/Ncol_${phenotype}.txt
    echo ${phenotype}
    echo "Sample size 1:${SAMPLE_SIZE1}" 
    echo "Sample size 2:${SAMPLE_SIZE2}"
    echo "Sample size 3:${SAMPLE_SIZE3}"
  
    paste "$IN"/${phenotype}_filtered.txt \
        <(awk -v n1="$SAMPLE_SIZE1" -v n2="$SAMPLE_SIZE2" -v n3="$SAMPLE_SIZE3" \
            'NR==1 {print "N"; next} NR>1 { print $2*n1 + $3*n2 + $4*n3 }' "$OUT"/${phenotype}_directions.txt) \
                > "$OUT"/Ncol_${phenotype}.txt
done

###########################################
# End
###########################################
