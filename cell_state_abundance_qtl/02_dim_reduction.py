# run this part as a separate script because it takes really long to run in jupyter notebook

import numpy as np
import pandas as pd
import multianndata as mad
import cna
import scanpy as sc

outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

madata = cna.read(f"{outdir}/data/scDataObject.h5ad")
# Define nearest neighbor graph and UMAP
print("Running PCA...")
sc.pp.pca(madata)
print("PCA finished!")

print("Generating nearest-neighbour graph...")
sc.pp.neighbors(madata)
print("Neighbourhood graph generation finished!")
print("Running UMAP ...")
sc.tl.umap(madata)
print("UMAP finished!")

# Define neighborhood abundance matrix (NAM) and NAM-PCs
print("Defining neighbourhood abundace matrix (NAM) ...")
cna.tl.nam(madata)
print("NAM finished!")
print("Saving outputs...")
madata.write(f"{outdir}/data/scDataObject.dimreduc.h5ad")
