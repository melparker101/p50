import scanpy as sc

# read in data
dataA = sc.read_h5ad('local.h5ad')

# extract raw count data
raw = dataA.raw.to_adata()
raw.to_df().iloc[:]
raw = raw.to_df()

# Test this
# Code for doing the python equivalent of R's table() apparently
# Don't really need this if using scanpy as it gives a summary at the bottom
df = pd.read_clipboard()
df.Gender.value_counts()   
