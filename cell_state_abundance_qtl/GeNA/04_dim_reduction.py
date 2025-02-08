import sys
import numpy as np
import pandas as pd
import multianndata as mad
import cna
import scanpy as sc
import scanpy.external as sce

celltype = sys.argv[1]
resolution = sys.argv[2]

outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

# ----
# Run PCA with Harmony to remove batch effects from the PCA components
# ----

madata = cna.read(f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.h5ad")
print("Running PCA...")
sc.pp.pca(madata)
print("PCA finished!")

print("Running Harmony...")
sce.pp.harmony_integrate(madata, "sequencing_library")
madata.obsm["X_pca"] = madata.obsm["X_pca_harmony"]
print("Harmony finished!")

# ----
# Define nearest neighbor graph and UMAP
# ----

print("Generating nearest-neighbour graph...")
sc.pp.neighbors(madata)
print("Neighbourhood graph generation finished!")

print("Running UMAP ...")
sc.tl.umap(madata)
print("UMAP finished!")

# ----
# Define neighborhood abundance matrix (NAM) and NAM-PCs
# ----

print("Defining neighbourhood abundace matrix (NAM) ...")
cna.tl.nam(
    madata
)  # note - run by GeNA as well (if not already present in the object) so can remove it here later
print("NAM finished!")
print("Saving outputs...")
madata.write(f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.h5ad")
