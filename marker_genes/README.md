### There are three Seurat functions we can use to find markers with (text below is from [Biostars](https://www.biostars.org/p/409790/); also see [hbctraining](https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html))

1. **FindMarkers** will find markers between two different identity groups - you have to specify both identity groups. This is useful for comparing the differences between two specific groups.

2. **FindAllMarkers** will find markers differentially expressed in each identity group by comparing it to all of the others - you don't have to manually define anything. Note that markers may bleed over between closely-related groups - they are not forced to be specific to only one group. This is what most people use (and likely what you want).

3. **FindConservedMarkers** will find markers that are conserved between two groups - this can be useful if you want to find markers that are conserved between a treated and untreated condition for a specific cell type or group of cells. It means they are differentially expressed compared to other groups, but have similar expression between the two groups you're actually comparing.

### Code for plots
https://lawrenson-lab.github.io/AtlasEndometriosis/figure2.html

### Order by pvalue and get top 10% of genes
- only.pos Only return positive markers (FALSE by default)

https://satijalab.org/seurat/reference/findmarkers 

# Info on datasets
### GSE118127
The object is filtered and preprocessed already. The data was normalised using the "LogNormalize" method, with some C++ code: https://github.com/johnmous/singleCell/blob/master/workflow.Rmd. 

### GSE202601
We only want to use data from the young patients, so this dataset needs filtering. It is also provided as raw counts, so we need to normalise it after filtering. We can use the Seurat function NormalizeData() with the "LogNormalize" and scale.factor = 10000 settings. Their original preprocessing code: https://github.com/satijalab/seurat/issues/678.

### GSE213216
The RDS file provided is already processed. It was normalised using the Seurat SCTransform() function. We want to filter to only include samples from patients 14 and 15 (young, European, and normal ovary samples). 
1) Filtering the processed (SCTransform) dataset
2) Extract the raw counts, filter, then normalise using "LogNormalize"
