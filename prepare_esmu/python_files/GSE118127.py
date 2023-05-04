############################################################
# This script is for running CELLEX on dataset GSE118127.
# melodyjparker14@gmail.com - Mar 2023
# Use the cellex conda environment.
############################################################

import scanpy as sc
import pandas as pd
import cellex

dataset_acc = 'GSE118127'

# Define data paths
dirIn = 'data/counts/' + dataset_acc + '/'  # 'counts/GSE118127/'
dirOut = 'data/esmu'

# Input files
input_file = 'local.h5ad'

# Output files
prefixData_sym = dataset_acc + '_sym'  # 'GSE118127_sym'
prefixData_ens = dataset_acc + '_ens'  # 'GSE118127_ens'

# Read in h5ad file to an annData object
adata = sc.read_h5ad(dirIn + input_file)

# Change any spaces in cell type names to underscores
adata.obs["cell_type"] = adata.obs["cell_type"].str.replace("[ ]","_")

# Extract metadata
metadata = adata.obs["cell_type"]
# metadata = pd.DataFrame(metadata)  # Not necessary

# Extract raw count data
raw = adata.raw.to_adata()
raw = raw.transpose()
mat = raw.to_df()
mat

# Check how many cells there are of each cell type
print(metadata.groupby('cell_type').cell_type.count())

### Run CELLEX ###
# Compute ESO object
eso = cellex.ESObject(data=mat, annotation=metadata, verbose=True)

# Compute ESw, ESw* and ESmu
eso.compute(verbose=True)

# Save and inspect results
# Save expression specificity mu and sd matrix for gene symbols
eso.save_as_csv(file_prefix=prefixData_sym, path=dirOut, verbose=True)
eso.results["esmu"].head()

# Map symbols to ensembl ids
cellex.utils.mapping.human_symbol_to_human_ens(eso.results["esmu"], drop_unmapped=True, verbose=True)
eso.results["esmu"].head()

# Save expression specificity matrix for gene ensembl ids
eso.save_as_csv(file_prefix=prefixData_ens, path=dirOut, verbose=True)

# Delete object to release memory
del eso
