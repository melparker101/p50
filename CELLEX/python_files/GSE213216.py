#######################################################
# Running CELLEX on dataset GSE213216
# melodyjparker14@gmail.com
# Use cellex conda env
# $PWD = p50
#######################################################

# change active cluster to cell_type in R

import scanpy as sc
import pandas as pd
import cellex

dataset_acc = 'GSE213216'

# Define data paths
dirIn = 'counts/' + dataset_acc + '/'
dirOut = 'esmu'

# Input files
input_file = 'counts_' + dataset_acc + '.h5ad'

# Output files
prefixData_sym = dataset_acc + '_sym'
prefixData_ens = dataset_acc + '_ens'

# Read in h5ad file to an annData object
adata = sc.read_h5ad(dirIn + input_file)

# Extract metadata
metadata = adata.obs["active.cluster"]
metadata = pd.DataFrame(metadata)

# Extract raw count data
raw = adata.raw.to_adata()
# raw = raw.transpose()
mat = raw.to_df()
mat

# NOT NORMALIESD

# Check how many cells there are of each cell type
print(metadata.groupby('active.cluster').active.cluster.count())

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
