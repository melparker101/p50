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

- sort out the sample size table - add into a separate script and make sure the table is written into premunge
