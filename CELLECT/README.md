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
3. Add rsid column to all sumstats using map file
4. Filter for MAF >1% in R
