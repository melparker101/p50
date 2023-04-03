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

```
