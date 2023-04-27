# Running CELLECT on ovary datasets
1. Munge sumstats
2. Clone the CELLECT github repository
- Clone the CELLECT git
``` bash
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
```
---
1. Create a map file (dbSNP/MarkerName_map_GRCh37.txt)
2. Add MarkerName column to hormones sumstats files (infertility already contains this)
3. Filter for MAF >1% in R
4. Add rsid column to all sumstats using map file

- write a script which changes yml files...
- maybe use multiple count datasets and sumstats in one file
- get sample size for infertility
- create 

``` bash
# pwd=CELLECT
conda activate <env_with_snakemake>
snakemake --use-conda -j -s cellect-ldsc.snakefile --configfile config.yml
```

### Things to note
1. Run CELLECT-magma and CELLECT-ldsc on rescomp the first time as it creates conda environments (internet connection is required for necessary package installation). After this, send slurm scripts off.
2. For CELLECT-ldsc, make sure the ldsc folder is properly downloaded - 'git clone' does not download this directory. Read the CELLECT github for more info.
3. For CELLECT-magma, use the 'keep p-val' option for munging - the summary stats must have a p-val column
