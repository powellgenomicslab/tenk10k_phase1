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

print(f"Creating MultiAnnData object for {celltype}...")

# constants
outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

# read in the latest tenk cohort
adata = sc.read(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_concatenated_harmony_filtered_min1000genes.h5ad",
    cache=True,
)

# filter only target cell type
adata = adata[adata.obs["wg2_scpred_prediction"] == celltype]

print("Cell type-specific adata:")
print(adata)

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

cell_meta.to_csv(f"{outdir}/data/metadata/{celltype}_cell_meta.csv")

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

sample_meta = (
    sample_meta.set_index("sample_id", drop=False)
    .loc[
        :, ~sample_meta.columns.str.startswith("sample_perm")
    ]  # remove the permutation columns
    .drop("sample_id", axis=1)
)

print("Sample metadata:")
print(sample_meta.head(3))
# for missing sex, use 0
sample_meta["sex"] = sample_meta["sex"].fillna(0)

madata = mad.MultiAnnData(X=cells_x_genes, obs=cell_meta, sampleid="id")
# Add all covariate information to d.samplem
madata.samplem = madata.samplem.join(
    sample_meta
)  # .join(batch_mapping) # don't add in the batch mappings for now
# check that you can actually save it with current metadata columns

print(madata)
madata.write(f"{outdir}/data/h5/{celltype}_scDataObject.h5ad")
