# Running CELLECT on ovary datasets
1. Munge sumstats
2. Clone the CELLECT github repository
- Clone the CELLECT git
``` bash
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
```
---
1. Add N column to sumstats ([add_N_col.sh](https://github.com/melparker101/p50/blob/main/CELLECT/add_N_col.sh))
2. Create a map file for chr:pos:a1_a2 to rsid (dbSNP/MarkerName_map_GRCh37.txt) ([create_snp_map.md](https://github.com/melparker101/p50/blob/main/CELLECT/create_snp_map.md) and [manipulate_map_file.sh](https://github.com/melparker101/p50/blob/main/CELLECT/manipulate_map_file.sh))
2. Add MarkerName column to hormones sumstats files (infertility sumstats already contains this) ([add_MarkerName_col_hormones.sh](https://github.com/melparker101/p50/blob/main/CELLECT/add_MarkerName_col_hormones.sh))
3. Filter for MAF >1% in R ([MAF_filter_all_sumstats.R](https://github.com/melparker101/p50/blob/main/CELLECT/MAF_filter_all_sumstats.R))
4. Add rsid column to all sumstats using map file
5. Munge using ldsc munge script

``` bash
# pwd=CELLECT
conda activate <env_with_snakemake>
snakemake --use-conda -j -s cellect-<method>.snakefile --configfile config.yml
```

### Things to note
1. Run CELLECT-magma and CELLECT-ldsc on rescomp the first time as it creates conda environments (internet connection is required for necessary package installation). After this, send slurm scripts off.
2. CELLECT creates these conda environments in a directory called '.snakemake'.
3. For CELLECT-ldsc, make sure the ldsc folder is properly downloaded - 'git clone' does not download this directory. Read the CELLECT github for more info.
4. For CELLECT-magma, use the 'keep p-val' option for munging - the summary stats must have a p-val column
5. An issue I came across was that CELLEX outputted a column of NAN for a certain cell type. This caused issues downstream when trying to run CELLECT. I am not sure what the solution is (temporary fix was to filter out that cell type from the data before running CELLEX, though we are losing information on that cell type)
6. When running CELLECT-magma, I was given an error message "ModuleNotFoundError: No module named 'statsmodels'; import statsmodels.api as sm". Although CELLECT is meant to install all packages automatically on the first run (apart from snakemake), I had to install this manually (conda install statsmodels) in the snakemake conda environment before attempting to run CELLECT-magma again - this fixed the issue.
