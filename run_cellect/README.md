# Run CELLECT

Clone the github repository. **Make sure --recurse-submodules argument is used to clone the ldsc directory**
``` bash
git clone --recurse-submodules https://github.com/perslab/CELLECT.git
```
Inside the CELLECT directory make a subdirectory **p50**; inside p50 create a subdirectory **data**. 

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


