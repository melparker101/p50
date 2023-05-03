########################################
# Make a text file with a marker column and a seperate direction column for each study
# 1 = SNP present in study, 0 = SNP not present in study
# melodyjparker14@gmail.com - Apr 23
########################################
# Run from p50 for now
# Use awk to edit as it is efficient

# Path for sumstats files
IN=data/sumstats/original
OUT=data/sumstats/other

########################################
# Create directions files
########################################

# Subset the sumstats data to only include the 'ID' and 'Direction' columns
for f in "$IN"/*EUR_filtered.txt
  do awk '{print $1, $12}' "$f" > "$OUT/$(basename -- "${f%_filtered.txt}")_directions.txt"
done

# The next part needs to be done manually as each file contains different studies

########################################
# LH 
########################################

# cohort  sample_size
# UKBB    10155
# EstBB   3353
# ALSPAC  3102

FILE_NAME=LH_F_EUR_directions.txt

"$FILE_NAME"

# Split up the characters of the 'Direction' column and name columns by study name
# Overwrite file
awk 'NR==1 {print $1, "UKBB", "EstBB", "ALSPAC"; next} NR==FNR {printf "%s ", $1; gsub(/.{1}/,"& ",$2); print $2}' "$OUT"/"$FILE_NAME" > tmp && mv tmp "$OUT"/"$FILE_NAME"

# Change ? to 0 and +/- to 1
awk '{gsub(/\?/,"0"); gsub(/[+-]/,"1"); print}' "$OUT"/"$FILE_NAME" > tmp && mv tmp "$OUT"/"$FILE_NAME"

########################################
# Oestradiol
########################################

# cohort  sample_size
# UKBB    49175
# EstBB   2740
# Pott    2607

FILE_NAME=Oestradiol_F_EUR_directions.txt

# Split up the characters of the 'Direction' column and name columns by study name
# Overwrite file
awk 'NR==1 {print $1, "UKBB", "EstBB", "Pott"; next} NR==FNR {printf "%s ", $1; gsub(/.{1}/,"& ",$2); print $2}' "$OUT"/"$FILE_NAME" > tmp && mv tmp "$OUT"/"$FILE_NAME"

# Change ? to 0 and +/- to 1
awk '{gsub(/\?/,"0"); gsub(/[+-]/,"1"); print}' "$OUT"/"$FILE_NAME" > tmp && mv tmp "$OUT"/"$FILE_NAME"

########################################
# Progesterone
########################################

# cohort  sample_size
# UKBB    1225
# EstBB   2683
# Pott    1261

FILE_NAME=Progesterone_F_EUR_directions.txt

# Split up the characters of the 'Direction' column and name columns by study name
# Overwrite file
awk 'NR==1 {print $1, "UKBB", "EstBB", "Pott"; next} NR==FNR {printf "%s ", $1; gsub(/.{1}/,"& ",$2); print $2}' "$OUT"/"$FILE_NAME" > tmp && mv tmp "$OUT"/"$FILE_NAME"

# Change ? to 0 and +/- to 1
awk '{gsub(/\?/,"0"); gsub(/[+-]/,"1"); print}' "$OUT"/"$FILE_NAME" > tmp && mv tmp "$OUT"/"$FILE_NAME"

########################################
# Testosterone
########################################

# cohort  sample_size
# UKBB    188196
# EstBB   3306
# ALSPAC  2380

FILE_NAME=Testosterone_F_EUR_directions.txt

# Split up the characters of the 'Direction' column and name columns by study name
# Overwrite file
awk 'NR==1 {print $1, "UKBB", "EstBB", "ALSPAC"; next} NR==FNR {printf "%s ", $1; gsub(/.{1}/,"& ",$2); print $2}' "$OUT"/"$FILE_NAME" > tmp && mv tmp "$OUT"/"$FILE_NAME"

# Change ? to 0 and +/- to 1
awk '{gsub(/\?/,"0"); gsub(/[+-]/,"1"); print}' "$OUT"/"$FILE_NAME" > tmp && mv tmp "$OUT"/"$FILE_NAME"

########################################
# FSH
########################################

# cohort  sample_size
# UKBB    11559
# EstBB   4457
# ALSPAC  3101

FILE_NAME=FSH_F_EUR_directions.txt

#########################################

# Define shell variables
FILE_NAME=FSH_F_EUR_directions.txt
COHORT1=UKBB
COHORT2=EstBB
COHORT3=ALSPAC

# Split up the characters of the 'Direction' column and name columns by study name
# Overwrite file
# Use -v to pass shell variables into awk command
awk -v c1="$COHORT1" -v c2="$COHORT2" -v c3="$COHORT3" \
  'NR==1 {print $1, c1, c2, c3; next} NR==FNR {printf "%s ", $1; gsub(/.{1}/,"& ",$2); print $2}' "$OUT"/"$FILENAME" \
    > tmp && mv tmp "$OUT"/"$FILENAME"

# Change ? to 0 and +/- to 1
awk '{gsub(/\?/,"0"); gsub(/[+-]/,"1"); print}' "$OUT"/"$FILENAME" > tmp && mv tmp "$OUT"/"$FILENAME"

#####################################
# Try this
FILE_NAME=FSH_F_EUR_directions.txt
COHORT1=UKBB
COHORT2=EstBB
COHORT3=ALSPAC

awk -v c1="$COHORT1" -v c2="$COHORT2" -v c3="$COHORT3" '{print $1, c1, c2, c3}' $OUT/$FILE_NAME | head -2

######################################


########################################
# End 
########################################
