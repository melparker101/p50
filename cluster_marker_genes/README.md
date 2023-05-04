# Find cluster marker genes
The following scripts find marker genes for the cell-type clusters of the single cell ovary datasets. Only clusters defined by the original study authors were used for each dataset.

### GSE118127
The object is filtered and preprocessed already. The data was normalised using the "LogNormalize" method, with some [C++ code](https://github.com/johnmous/singleCell/blob/master/workflow.Rmd).

### GSE202601
We only want to use data from the young patients, so this dataset needs filtering. It is also provided as raw counts, so we need to normalise it after filtering. We can use the Seurat function NormalizeData() with the "LogNormalize" and scale.factor = 10000 settings. Their original preprocessing code: https://github.com/satijalab/seurat/issues/678.

### GSE213216
The RDS file provided is already processed. It was normalised using the Seurat SCTransform() function. We want to filter to only include samples from patients 14 and 15 (young, European, and normal ovary samples). 
1) Filtering the processed (SCTransform) dataset
2) Extract the raw counts, filter, then normalise using "LogNormalize"

### URLs
https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html
https://satijalab.org/seurat/reference/findmarkers 

