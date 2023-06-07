# This script reads in a count matrix and metadata, runs cellex, and outputs esmu and essd
# Use the cellex conda env

import scanpy as sc
import pandas as pd
import cellex

# Define data paths
dirIn = 'data/counts/GSE202601/'
dirOut='data/esmu'

# Input files
matrix_file = 'counts_GSE202601.h5ad'
metadata_file = 'metadata_GSE202601.csv'

# Output files
prefixData_sym='GSE202601_sym'
prefixData_ens='GSE202601_ens'

# Read in data
adata = sc.read_h5ad(dirIn + matrix_file)
adata = adata.transpose()
mat = adata.to_df()
mat

metadata = pd.read_csv(dirIn + metadata_file, index_col=0)
metadata

# Check how many cells there are of each cell type
print(metadata.groupby('cell_type').cell_type.count())

# Run CELLEX
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

# Other save options
# eso.save_as_csv(keys=["all"], verbose=True)  # Save all results
# eso.results["esmu"].to_csv("esmu/" + "GSE202601.esmu.csv.gz")

# Delete object to release memory
del eso
