import sys
import numpy as np
import pandas as pd
import multianndata as mad
import cna
import scanpy as sc
import scanpy.external as sce
from matplotlib import pyplot as plt


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
# NOTE: might be able to tune this better by using fewer PC's in the NN graph. GeNA paper used 20 PCs, default is 50 PCs
sc.pp.pca(madata)
print("PCA finished!")

print("Running Harmony...")
# NOTE: Maybe try with harmony sigma parameter 0.2 to encourange softer clustering for the major cell types
sce.pp.harmony_integrate(madata, "sequencing_library")
madata.obsm["X_pca"] = madata.obsm["X_pca_harmony"]
print("Harmony finished!")

# ----
# Define nearest neighbor graph and UMAP
# ----

print("Generating nearest-neighbour graph...")
sc.pp.neighbors(madata)
print("Neighbourhood graph generation finished!")

if celltype != "ALL":
    print("Running UMAP ...")
    sc.tl.umap(madata)
    print("UMAP finished!")

# ----
# Define neighborhood abundance matrix (NAM) and NAM-PCs
# ----

print("Defining neighbourhood abundace matrix (NAM) ...")
cna.tl.nam(
    madata
)  # note - run by GeNA as well (if not already present in the object) so can remove it here later if need to speed up
print("NAM finished!")
print("Saving outputs...")
madata.write(f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.h5ad")

# plot umaps
color_vars = ["wg2_scpred_prediction", "cohort"]
with plt.rc_context({"figure.figsize": (5, 5)}):
    sc.pl.umap(madata, color=color_vars, frameon=False, ncols=2)
plt.savefig(f"{outdir}/figures/{resolution}/{celltype}_umap.png")
