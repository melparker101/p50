# Run from p50 for now

# Path for sumstats files
ORIG_SUMSTATS=data/sumstats/original

# Subset the sumstats data to only include the 'MarkerName' and 'Direction' columns
awk '{print $1, $11}' "$ORIG_SUMSTATS"/female_infertility_analysis1_UKBB_Finngen_EstBB_noMACfilter_March20231.out > "$ORIG_SUMSTATS"/Infertility1_F_EUR_directions.txt

# Split up the characters of the 'Direction' column and overwrite file
awk -F" " '{printf "%s\t", $1; gsub(/.{1}/,"& ",$2); print $2}' "$ORIG_SUMSTATS"/Infertility1_F_EUR_directions.txt > tmp && mv tmp "$ORIG_SUMSTATS"/Infertility1_F_EUR_directions.txt

# Infertility1_F_EUR.txt
