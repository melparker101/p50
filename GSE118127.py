import scanpy as sc

# read in data
dataA = sc.read_h5ad('local.h5ad')

# extract raw count data
raw = dataA.raw.to_adata()
raw.to_df().iloc[:]

# Test this
# Code for doing the python equivalent of R's table() apparently
df = pd.read_clipboard()
df.Gender.value_counts()   
