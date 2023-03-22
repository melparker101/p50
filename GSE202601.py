import loompy # needed for importing data for this tutorial
import numpy as np # needed for formatting data for this tutorial
import pandas as pd # needed for formatting data for this tutorial
import cellex

dirOut = "cellex_demo_out" # output directory for results and plots
prefixData = "mousebrain_vascular_cells" # prefix to prepend to files

pathData = "GSE202601.loom"
nameAnno = "ClusterName" # metadata annotation column attribute name
nameId = "CellID" # metadata cell id column attribute name
nameClass = "Class"

with loompy.connect(pathData) as ds:
    rows = (ds.row_attrs["Accession"])
    cols = (ds.col_attrs[nameId])
    data = pd.DataFrame(ds[:, :], index=rows, columns=cols)
    n_cells_total = data.shape[1]
    
    
pathData = "GSE202601.loom"

ds = loompy.connect("GSE202601.loom")
# do something with the connection object ds
ds.close()
    
  # not working... we could read in with scanpy then extract count matrix?
    
 or just read in as a normal h5 file using h5py
