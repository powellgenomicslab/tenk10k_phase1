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

# celltype = "Other"
# resolution = "major_cell_types"

# constants
outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

print(f"Preprocessing anndata for {celltype}...")

# read in the latest tenk cohort
# had to split previous steps to different script due to memory issues
adata = sc.read(
    f"{outdir}/data/h5/{resolution}/{celltype}_scanpy.h5ad",
    cache=True,
)

if celltype == "ALL":
    major_cell_type_mapping = {
        "B_intermediate": "B",
        "B_memory": "B",
        "B_naive": "B",
        "Plasmablast": "B",
        "NK": "NK",
        "NK_CD56bright": "NK",
        "NK_Proliferating": "NK",
        "CD8_Naive": "CD8_T",
        "CD8_Proliferating": "CD8_T",
        "CD8_TCM": "CD8_T",
        "CD8_TEM": "CD8_T",
        "CD4_CTL": "CD4_T",
        "CD4_Naive": "CD4_T",
        "CD4_Proliferating": "CD4_T",
        "CD4_TCM": "CD4_T",
        "CD4_TEM": "CD4_T",
        "Treg": "CD4_T",
        "dnT": "Unconventional_T",
        "gdT": "Unconventional_T",
        "ILC": "Unconventional_T",
        "MAIT": "Unconventional_T",
        "pDC": "Dendritic",
        "cDC1": "Dendritic",
        "cDC2": "Dendritic",
        "ASDC": "Dendritic",
        "CD14_Mono": "Monocyte",
        "CD16_Mono": "Monocyte",
        "HSPC": "Other",
    }

    # add major cell type to single cell metadata

    adata.obs["major_cell_type"] = [
        major_cell_type_mapping[ct] for ct in adata.obs["wg2_scpred_prediction"]
    ]


# adata.X = adata.layers["counts"].copy()  # so that it doesn't error on a re-run
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

print("Cell type-specific adata:")
print(adata)

# overwrite h5ad
# rerun previous script before rerunning above steps otherwise it will error
adata.write(f"{outdir}/data/h5/{resolution}/{celltype}_scanpy.h5ad")

# ----------------------
# get the cell metadata
# ----------------------

cell_meta = (
    adata.obs[["cpg_id", "sequencing_library", "cohort", "wg2_scpred_prediction"]]
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
