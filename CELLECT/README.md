## Using CELLEX and CELLECT on single cell RNA-seq ovary datasets with infertility GWAS summary statistics
To use CELLEX and CELLECT, follow the instructions on their github repositories. Once the CELLECT directory is cloned from their github, create a subdirectory **p50** for this project.

### Directory structure
The basic directory structure for the p50 infertility project work:
```
p50
|-- CELLECT_OUT_p50
|   |-- CELLECT-GENES
|   |-- CELLECT-LDSC
|   `-- CELLECT-MAGMA
|-- cluster_markers
|   |-- GSE118127
|   |-- GSE189960
|   |-- GSE192722
|   |-- GSE202601
|   |-- GSE206143
|   `-- GSE213216
|-- data
|   |-- counts
|   |-- esmu
|   `-- sumstats
|-- dbSNP
|-- logs
`-- plots
```
- **CELLECT_OUT_p50** - Created when CELLECT is run. Contains CELLECT output files.
- **cluster_markers** - Store cluster_marker_genes output files here.
- **data** - Store input data for CELLEX and CELLECT here.
- **dbSNP** - Store the MarkerName to RSID map file here.
- **logs** - For logs.
- **plots** - Store plots generated from CELLECT results here. 

### Pipeline
Once the directory structure is set up, follow the pipeline below.
1. **Download data**  
We need scRNA-seq count data and the corresponding cell type annotations metadata. We also need GWAS summary statistics (in-house). See [datasets](https://github.com/melparker101/p50/tree/main/datasets) for more information.  
2. **Set up environments**  
Download required packages/create the recommended conda environments. More information is given in [set_up](https://github.com/melparker101/p50/tree/main/set_up).
3. **Prepare ESMU files (run CELLEX)**  
Using the counts and cell type annotations metadata as input, we use [CELLEX](https://github.com/perslab/CELLEX) to produce expression specificity files (ESMU). See [prepare_esmu](https://github.com/melparker101/p50/tree/main/prepare_esmu) for R and python code used to prepare data and run CELLEX.
4. **Prepare sumstats file**  
Use the pipeline provided in [prepare_sumstats](https://github.com/melparker101/p50/tree/main/prepare_sumstats) to prepare the GWAS summary statistics for input to CELLECT.  
5. **Run CELLECT**  
Using munged summary stats and ESMU files as input, we use [CELLECT](https://github.com/perslab/CELLECT/wiki/CELLECT-LDSC-Tutorial) to prioritise etilogical cell types. Use the **config_p50.yml** file provided. See [run_cellect](https://github.com/melparker101/p50/tree/main/run_cellect). 
6. **Visualisation**  
Use R to visualise the results. See [visualisation](https://github.com/melparker101/p50/tree/main/visualisation).  



source ~/.bashrc so that conda can be activated. make sure the conda path is in bashrc.
