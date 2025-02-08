import sys
import numpy as np
import pandas as pd
import multianndata as mad
import cna
import scanpy as sc
import matplotlib.pyplot as plt

np.random.seed(0)

# get cell type from cmdline arg 1
celltype = sys.argv[1]
resolution = sys.argv[2]

celltype = "CD4_T"
resolution = "major_cell_types"

# constants
outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

# celltype = "NK"
# resolution = "major_cell_types"

print(f"Preprocessing anndata for {celltype}...")

# read in the latest tenk cohort
adata = sc.read(
    f"{outdir}/data/h5/{resolution}/{celltype}_scanpy.h5ad",
    cache=True,
)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_ct = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

print("Cell type-specific adata:")
print(adata)

# overwrite h5ad
adata.write(f"{outdir}/data/h5/{resolution}/{celltype}_scanpy.h5ad")

# ----------------------
# get the cell metadata
# ----------------------

cell_meta = (
    adata.obs[["cpg_id", "sequencing_library"]]
    .set_index(adata.obs.index)
    .rename(columns={"cpg_id": "id"})
)

print("Cell metadata:")
print(cell_meta.head(3))

# ----------------------
# get the counts matrix
# ----------------------

cells_x_genes = adata.to_df()
# export counts to csv
print("Cells x genes matrix:")
print(cells_x_genes.shape)
print(cells_x_genes.head(3))

# ----------------------
# get the sample metadata (covariates)
# ----------------------

sample_meta = pd.read_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/saige-qtl_tenk10k-genome-2-3-eur_input_files_241210_covariates_sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv"
)
sample_meta = sample_meta[sample_meta["sample_id"].isin(adata.obs["cpg_id"])]

sample_meta = sample_meta.set_index("sample_id", drop=True)

print("Sample metadata:")
print(sample_meta.head(3))
# for missing sex, use 0
sample_meta["sex"] = sample_meta["sex"].fillna(0)

print(f"Creating MultiAnnData object for {celltype}...")

madata = mad.MultiAnnData(X=cells_x_genes, obs=cell_meta, sampleid="id")
# Add all covariate information to d.samplem
madata.samplem = madata.samplem.join(sample_meta)

print(madata)
madata.write(f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.h5ad")
