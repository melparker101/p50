########################################
# Make a text file with a marker column and a seperate direction column for each study
# melodyjparker14@gmail.com - Apr 23
########################################
# Run from p50 for now
# Use awk to edit as it is efficient

# Path for sumstats files
IN=data/sumstats/original
OUT=data/sumstats/not_using

# Subset the sumstats data to only include the 'MarkerName' and 'Direction' columns
awk '{print $1, $11}' "$IN"/female_infertility_analysis1_UKBB_Finngen_EstBB_noMACfilter_March20231.out > "$OUT"/Infertility1_F_EUR_directions.txt

# Split up the characters of the 'Direction' column and name columns by study name
# Overwrite file
awk 'NR==1 {print $1, "FinnGen", "UKBB", "EstBB"; next} NR==FNR {printf "%s ", $1; gsub(/.{1}/,"& ",$2); print $2}' "$OUT"/Infertility1_F_EUR_directions.txt > tmp && mv tmp "$OUT"/Infertility1_F_EUR_directions.txt

########################################
# End 
########################################
