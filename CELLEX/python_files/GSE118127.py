############################################################
# This script is for running CELLEX on dataset GSE118127.
# melodyjparker14@gmail.com
# Use the cellex conda environment.
############################################################

import scanpy as sc
import pandas as pd
import cellex

# Read in data
adata = sc.read_h5ad('local.h5ad')

# Extract metadata
metadata = adata.obs["cell_type"]
metadata = pd.DataFrame(metadata)

# Extract raw count data
raw = adata.raw.to_adata()
raw = raw.transpose()
mat = raw.to_df()
mat

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
