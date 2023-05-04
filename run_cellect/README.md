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



