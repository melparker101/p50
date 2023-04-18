### There are three Seurat functions we can use to find markers with (see [Biostars](https://www.biostars.org/p/409790/) and [hbctraining](https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_marker_identification.html))

1. FindMarkers will find markers between two different identity groups - you have to specify both identity groups. This is useful for comparing the differences between two specific groups.

2. FindAllMarkers will find markers differentially expressed in each identity group by comparing it to all of the others - you don't have to manually define anything. Note that markers may bleed over between closely-related groups - they are not forced to be specific to only one group. This is what most people use (and likely what you want).

3. FindConservedMarkers will find markers that are conserved between two groups - this can be useful if you want to find markers that are conserved between a treated and untreated condition for a specific cell type or group of cells. It means they are differentially expressed compared to other groups, but have similar expression between the two groups you're actually comparing.
