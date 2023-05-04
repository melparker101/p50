# Prepare sumstats
GWAS summary statistics need to be in a specific format to use as input for CELLECT. The following pipeline was used to manipulate the files.

### Pipeline
1. Add N column to sumstats ([add_N_col.sh](https://github.com/melparker101/p50/blob/main/CELLECT/add_N_col.sh))
2. Create a map file for chr:pos:a1_a2 to rsid (dbSNP/MarkerName_map_GRCh37.txt) ([create_snp_map.md](https://github.com/melparker101/p50/blob/main/CELLECT/create_snp_map.md) and [manipulate_map_file.sh](https://github.com/melparker101/p50/blob/main/CELLECT/manipulate_map_file.sh))
2. Add MarkerName column to hormones sumstats files (infertility sumstats already contains this) ([add_MarkerName_col_hormones.sh](https://github.com/melparker101/p50/blob/main/CELLECT/add_MarkerName_col_hormones.sh))
3. Filter for MAF >1% in R ([MAF_filter_all_sumstats.R](https://github.com/melparker101/p50/blob/main/CELLECT/MAF_filter_all_sumstats.R))
4. Add rsid column to all sumstats using map file ([add_rsid.sh](https://github.com/melparker101/p50/blob/main/CELLECT/add_rsid.sh))
5. Munge using ldsc munge script ([munge_sumstats.sh](https://github.com/melparker101/p50/blob/main/CELLECT/munge_sumstats.sh))

### Directory structure 
Structure of directories for sumstats after running the pipeline
``` text
p50/data/sumstats/
|-- cohorts
|-- munged
|-- original
`-- other
```
The **munged** directory will contain the final sumstats files which will be used as input for CELLECT.
