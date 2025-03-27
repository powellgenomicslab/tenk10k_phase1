import sys
import numpy as np
import pandas as pd
import cna
import multianndata as mad
import os

celltype = sys.argv[1]
resolution = sys.argv[2]

# resolution = "major_cell_types"
# analysis_name = "no_expr_pc_covars"
# celltype = "NK"

# input files
madata_file = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.pca.h5ad"
madata = cna.read(madata_file)

umap_df = pd.DataFrame(
    madata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=madata.obs.index
)

umap_df.to_csv(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/umaps/{resolution}/{celltype}_umap.csv"
)
