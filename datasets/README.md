# Datasets

### Ovary scRNA-seq datasets
The GSE .md files in this sub directory contain general information about each dataset used.
- GSE202601
- GSE118127 
- GSE213216 *
- GSE189960 
- GSE192722 
- GSE206143 

\* This dataset was not used for CELLECT! It does not contain cell types of interest and also did not work with CELLEX.

### GWAS Sumstats (build 37/hg19)
Move the in-house GWAS summary statistic files into p50/data/sumstats/original.
- female_infertility_analysis1_UKBB_Finngen_EstBB_noMACfilter_March20231.out
- FSH_F_EUR_filtered.txt
- Oestradiol_F_EUR_filtered.txt
- LH_F_EUR_filtered.txt          
- Progesterone_F_EUR_filtered.txt
- Testosterone_F_EUR_filtered.txt
- Testosterone_sex_comb_EUR_filtered.txt

### Directory structure
Layout of the data directory:
```
p50
`-- data
    |-- counts
        |-- GSE118127
        |-- GSE189960
        |-- GSE192722
        |-- GSE202601
        |-- GSE206143
        `-- GSE213216
    `-- sumstats
        `-- original
```
