# Find cluster marker genes
The following scripts find marker genes for the cell-type clusters of the single cell ovary datasets. Only clusters defined by the original study authors were used for each dataset.

| Acc No.   | Paper URL                                                                           | Dataset URL                                                   | Cells  | Number of clusters | Github                                                                           | Cell-type annotations provided? |
|-----------|-------------------------------------------------------------------------------------|---------------------------------------------------------------|--------|--------------------|----------------------------------------------------------------------------------|---------------------------------|
| GSE189960 | [ncbi](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8961562/)                                        | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189960)  | 14,592 | 13                 | -                                                                                | No                              |
| GSE202601 | [biorxiv](https://www.biorxiv.org/content/biorxiv/early/2022/05/19/2022.05.18.492547.full.pdf) | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202601)  | 42,568 | 8                  | https://github.com/ChenJin2020/The-regulatory-landscapes-of-human-ovarian-ageing | Yes                             |
| GSE192722 | [nature](https://www.nature.com/articles/s42003-022-04384-8)                                  | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192722)  | 48,147 | 22,6               | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192722                     | No                              |
| GSE118127 | [ncbi](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6639403/)                               | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118127)  | 20,676 | 19,5               | https://github.com/johnmous/singleCell                                           | Yes                             |
| GSE213216 | [nature](https://www.nature.com/articles/s41588-022-01254-1)                                  | [Cedars](https://cedars.app.box.com/s/1ks3eyzlpnjbrseefw3j4k7nx6p2ut02); [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213216) | 22,219 | 9                  | https://github.com/lawrenson-lab/AtlasEndometriosis                              | Yes                             |
| GSE206143 | [wiley](https://faseb.onlinelibrary.wiley.com/doi/10.1096/fj.202201746RR)                    | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206143)  | 7609   | 19,13              | https://github.com/nurungji82/scRNA-seq_of_IVF_samples                           | No                              |

### GSE118127
The object is filtered and preprocessed already. The data was normalised using the "LogNormalize" method, with some [C++ code](https://github.com/johnmous/singleCell/blob/master/workflow.Rmd).

### GSE202601
We only want to use data from the young patients, so this dataset needs filtering. It is also provided as raw counts, so we need to normalise it after filtering. We can use the Seurat function NormalizeData() with the "LogNormalize" and scale.factor = 10000 settings. Their original preprocessing code: https://github.com/satijalab/seurat/issues/678.

### GSE213216
The RDS file provided is already processed. It was normalised using the Seurat SCTransform() function. We want to filter to only include samples from patients 14 and 15 (young, European, and normal ovary samples). 
1) Filtering the processed (SCTransform) dataset
2) Extract the raw counts, filter, then normalise using "LogNormalize"

### GSE189960  
Cell-type metadata was not provided. Clustering and cell-type annotations were replicated as closely as possible, manually, following the methods provided in the paper and github code if availible. We do not know the ethnicity of the donors, so all samples were used.

### GSE192722
Cell-type metadata was not provided. Clustering and cell-type annotations were replicated as closely as possible, manually, following the methods provided in the paper and github code if availible. We do not know the ethnicity of the donors, so all samples were used.

### GSE206143
Cell-type metadata was not provided. Clustering and cell-type annotations were replicated as closely as possible, manually, following the methods provided in the paper and github code if availible. We do not know the ethnicity of the donors, so all samples were used.

### URLs
https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html
https://satijalab.org/seurat/reference/findmarkers 

