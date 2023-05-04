# Run CELLECT
CELLECT takes two inputs:
- Expression specificity CSV (ESMU). Use CELLEX to create this. Use ensembl ids.
- GWAS summary statistics. Use this pipeline to prepare the sumstats.

This is how the data should be organised to use the config_p50.yml file:
```
p50/data/
|-- counts
|   |-- GSE118127
|   |-- GSE202601
|   `-- GSE213216
|-- esmu
`-- sumstats
    |-- cohorts
    |-- munged
    |-- original
    `-- other
```

There are three versions of CELLECT:
- [CELLECT LDSC](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial)
- [CELLECT MAGMA](https://github.com/perslab/CELLECT/wiki/CELLECT-MAGMA-Tutorial)
- [CELLECT GENES](https://github.com/perslab/CELLECT/wiki/CELLECT-GENES-Tutorial)

### Things to note
1. Run CELLECT-magma and CELLECT-ldsc on rescomp the first time as it creates conda environments (internet connection is required for necessary package installation). After this, send slurm scripts off.
2. CELLECT creates these conda environments in a directory called '.snakemake'.
3. For CELLECT-ldsc, make sure the ldsc folder is properly downloaded - 'git clone' does not download this directory. Make sure --recurse-submodules argument is used to clone the ldsc directory.
``` bash
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
```
4. For CELLECT-magma, use the 'keep p-val' option for munging - the summary stats must have a p-val column
5. An issue I came across was that CELLEX outputted a column of NAN for a certain cell type. This caused issues downstream when trying to run CELLECT. I am not sure what the solution is (temporary fix was to filter out that cell type from the data before running CELLEX, though we are losing information on that cell type)
6. When running CELLECT-magma, I was given an error message "ModuleNotFoundError: No module named 'statsmodels'; import statsmodels.api as sm". Although CELLECT is meant to install all packages automatically on the first run (apart from snakemake), I had to install this manually (**conda install statsmodels**) in the snakemake conda environment before attempting to run CELLECT-magma again - this fixed the issue.
7. For CELLECT-genes, use this temporary bug fix: https://github.com/perslab/CELLECT/issues/81

### CELLECT output directory structure
```
p50/CELLECT_OUT_p50
|-- CELLECT-GENES
|   |-- logs
|   |-- out
|   `-- results
|-- CELLECT-LDSC
|   |-- logs
|   |-- out
|   |   `-- prioritization
|   |-- precomputation
|   |   |-- GSE118127
|   |   |   `-- per_annotation
|   |   |-- GSE202601
|   |   |   `-- per_annotation
|   |   |-- GSE213216
|   |   |   `-- per_annotation
|   |   |-- bed
|   |   `-- control.all_genes_in_dataset
|   `-- results
`-- CELLECT-MAGMA
    |-- logs
    |-- out
    |   `-- prioritization
    |-- precomputation
    |   |-- FSH_F_EUR
    |   |-- Infertility1_F_EUR
    |   |-- LH_F_EUR
    |   |-- Oestradiol_F_EUR
    |   |-- Progesterone_F_EUR
    |   |-- Testosterone_F_EUR
    |   `-- Testosterone_sex_comb_EUR
    `-- results
```
