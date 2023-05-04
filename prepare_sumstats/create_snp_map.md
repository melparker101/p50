# Create map file for GRCh37
These instructions are for creating a SNP map file for GRCh37, but the same process could be followed to create a map for GRCh38. First, download the required VCF file from website.

``` text
/VCF               RefSNP VCF files for GRC (Genome Reference Consortium) human assembly
                   37 (GCF_000001405.25) and 38 (GCF_000001405.40). Files are compressed
                   by bgzip and with the tabix index.
```

``` bash
# Load modules
module load BCFtools/1.17-GCC-12.2.0

# Download the desired VCF file.
mkdir dbSNP
cd dbSNP
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz

# Download the corresponding index file
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi

# Alternatively, index the VCF file manually (the VCF file must be zipped)
# bcftools index -t GCF_000001405.25.gz

# We don't need to unzip, but if we did use this
# bgzip -d GCF_000001405.25.gz

# Take a look at the relevant VCF columns
bcftools query -f '%CHROM %POS %ID %REF %ALT\n' GCF_000001405.25.gz | head -3
```

The chromosomes are in the wrong format.
``` text
NC_000001.10 10001 rs1570391677 T A,C
NC_000001.10 10002 rs1570391692 A C
NC_000001.10 10003 rs1570391694 A C
```
We need to convert these... Next, download the [GRCh37 assembly report](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt) 
to manipulate and use as a reference (See [biostars](https://www.biostars.org/p/410789/) for code).

``` bash
# Load modules
module load BCFtools/1.17-GCC-12.2.0

# Download assembly report
report_dir='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405'
wget -N "${report_dir}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt"

# Grab the useful columns
for k in *assembly_report.txt
  do
    out=$(echo $k | sed 's/.txt/.chrnames/')
    grep -e '^[^#]' $k | awk '{ print $7, $1 }' > $out
done

# Annotate (this step takes a while)
bcftools annotate \
--rename-chrs GCF_000001405.25_GRCh37.p13_assembly_report.chrnames \
GCF_000001405.25.gz \
-Oz -o corrected_GCF_000001405.25.gz

# Check that the conversion was successful
bcftools query -f '%CHROM %POS %ID %REF %ALT\n' corrected_GCF_000001405.25.gz | head -3

# Create a map file
bcftools query -f '%CHROM %POS %ID %REF %ALT\n' corrected_GCF_000001405.25.gz > GRCh37_GCF_000001405.25_map.txt
```

The map file should look like this (chr, pos, rsid, allele1, allele2):
``` text
1 10001 rs1570391677 T A,C
1 10002 rs1570391692 A C
1 10003 rs1570391694 A C
```
The file has 1,118,410,664 lines! Filter the allele columns to only keep lines where they both contain one character only.

``` bash
awk 'length($4)==1 && length($5)==1' GRCh37_GCF_000001405.25_map.txt > map_filtered.txt
```

```bash
# awk version of table() from R
awk -F" " '{arr[$4]++}END{for (a in arr) print a, arr[a]}' test2.txt

# Find all chrm numbers in map file
awk -F" " '{arr[$1]++}END{for (a in arr) print a, arr[a]}' map_sv_filtered.txt > table.txt

# there are patches included
sed '/H/d' map_sv_filtered.txt > map_filtered.txt

# test again
awk -F" " '{arr[$1]++}END{for (a in arr) print a, arr[a]}' map_filtered.txt > table_f.txt
cat table_f.txt
we are now left with chromosomes 1-22, X, Y and MT

# Do the same for LAuras sumstats
awk -F: '{arr[$1]++}END{for (a in arr) print a, arr[a]}' female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out > table_laura.txt

# Remove patches (any chrom starting with H), mitochondrial and the Y chomosome.
sed -e '/H/d' -e '/MT/d' -e '/Y/d' map_sv_filtered.txt > map_filtered2.txt

# check this
awk -F" " '{arr[$1]++}END{for (a in arr) print a, arr[a]}' map_filtered2.txt > table_f2.txt
wc -l map_filtered2.txt

# Do the same for samvidas
awk -F: '{arr[$1]++}END{for (a in arr) print a, arr[a]}' FSH_F_EUR_filtered.txt > table_samvida.txt
cat table_samvida.txt
# we have chr1-23

# remove doubles and M and Y

# Make a MarkerName map file - this will take a while so use a few cores or send a script off
awk '{if ($4 < $5) {a1=$4; a2=$5} else {a1=$5; a2=$4}; print $1":"$2":"a1"_"a2,$3}' map_filtered2.txt > MarkerName_map_GRCh37.txt

# Separate into chr files
awk -F: '{if ($1==1)print $0}' MarkerName_map_GRCh37.txt > map_by_chrom.37/map_chr1.txt
awk
awk -F\| '{print>$1}' file1

# Add to separate files per chrom (>> isn't necessary as it is keeps the files open)
awk -F: '{print > "maps/map_chr"$1".txt"}' MarkerName_map_GRCh37.txt

# first, sort the sumstats in numerical order
sort -n female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out > female_infertility_UKBB.txt

awk -F'\t' 'NR==FNR{a[$1];next} $2 in a{print $2, $1}' test_ss.txt ../dbSNP/MarkerName_map_GRCh37.txt

awk -F'\t' 'NR==FNR{a[$1];next} ($2 in a) && !seen[$2]++{print $2, $16}' Accessions gene2accession

awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1]=$0; next} $1 in a {print a[$1],$0}' test_ss.txt ../dbSNP/MarkerName_map_GRCh37.txt > test_joined.txt

# Might be easier to do it using awk
awk -F'\t' 'NR==FNR{a[$1];next} ($2 in a) && !seen[$2]++{print $2, $16}' MarkerName_map_GRCh37.txt female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out
```

```R

library(stringr)
word(ss$MarkerName,1,sep = ":")

while(length(sumstats_row <- read.table(sumstats_file, header=FALSE, sep='\t', nrows=1)) > 0) {
    # Extract the rows from mapping_reader that match the current chromosome
    mapping_rows <- subset(mapping_reader, word(mapping_reader[1],1,sep = ":") == sumstats_row[1])

    # Loop through each row in the mapping file until we find a match
    rsid <- NULL
    for(i in 1:nrow(mapping_rows)) {
        mapping_row <- mapping_rows[i,]
        if (mapping_row[2] == sumstats_row[2]) {
            rsid <- mapping_row[3]
            break
        }
    }
```

```
# Open the mapping file
mapping_file <- file('MarkerName_map_GRCh37.txt', 'r')
mapping_reader <- read.table(mapping_file, header=FALSE, sep='\t')
############################

ss_string <- '../sumstats/female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out'

# Open the sumstats file
sumstats_file <- file('../sumstats/female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out', 'r')
sumstats_reader <- read.table(sumstats_file, header=TRUE, sep='\t')

# fread is much much faster - read in sumstats with this

# Open the output file
output_file <- file('merged.txt', 'w')
output_writer <- write.table(x=NULL, file=output_file, sep='\t', quote=FALSE, row.names=FALSE)

# Write the header row to the output file
header_row <- names(sumstats_reader)
header_row[length(header_row) + 1] <- 'rsid'
write.table(x=header_row, file=output_file, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  
while(length(sumstats_row <- read.table(sumstats, header=FALSE, sep='\t', nrows=1)) > 0) {
    # Extract the rows from mapping_reader that match the current chromosome
    chr <- word(sumstats$MarkerName, 1, sep = ":")
    map <- paste0('maps/map_chr",chr,".txt)
    
    mapping_file <- file(map, 'r')
    mapping_reader <- read.table(mapping_file, header=FALSE, sep='\t')

    # Loop through each row in the mapping file until we find a match
    rsid <- NULL
    for(i in 1:nrow(mapping_reader)) {
        mapping_row <- mapping_rows[i,]
        if (mapping_row[1] == sumstats_row[1]) {
            rsid <- mapping_row[2]
            break
        }
    }
       
    # Write the merged row to the output file
    sumstats_row[length(sumstats_row) + 1] <- rsid
    write.table(x=sumstats_row, file=output_file, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
}
    }
  
  
  
  ##################
  
  # Loop through each row in the sumstats file
for (i in 1:nrow(sumstats)) {
    # Subset the mapping file to find a match
    chr <- word(sumstats$MarkerName, 1, sep = ":")
    map <- paste0('maps/map_chr",chr,".txt)
    mapping <- fread(sprintf("awk '$1' ../sumstats/female_infertility_analysis1_UKBB_Finngen_EstBB_GandH_noMACfilter_March20231.out", sumstats$chr[i], sumstats$pos[i]))
    if (nrow(mapping) > 0) {
        rsid <- mapping$rsid[1]
    } else {
        rsid <- NA
    }

    # Write the merged row to the output file
    merged_row <- c(sumstats[i,], rsid)
    write.table(x=merged_row, file=output_file, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
}
  
  
# Close the files
close(mapping_file)
close(sumstats_file)
close(output_file)

60


closeAllConnections()

```
