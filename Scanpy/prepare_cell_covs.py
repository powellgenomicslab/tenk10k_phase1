import os
import sys
import re
import pandas as pd
import scanpy as sc
import scanpy.external as sce

celltype = sys.argv[1]

# Specify which directory the files generated here will be saved to
output_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/cpg_cell_covs_filtered"

if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# Specify which files this script will generate
# CSV for expression covariates (principal components after Harmony correction)
pcs_out_file = f"{output_dir}/{celltype}_expression_pcs_no_harmony_with_regression.csv"
# if os.path.exists(pcs_out_file):
#   sys.exit("File already exists!")

# Load combined AnnData object
input_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries"

# NOTE: actually run with:
# input_file = f"{input_dir}/300_libraries_concatenated_filtered.h5ad"
# This is a copy of: 298_libraries_cell_metadata_filtered_min1000genes.csv (updated file name)
# confirmed same number of cells
input_file = f"{input_dir}298_libraries_cell_metadata_filtered_min1000genes.csv"
adata = sc.read(input_file)

# Extract cell type specific expression
adata_ct = adata[adata.obs["wg2_scpred_prediction"] == celltype]

# normalisation (per cell type)
sc.pp.normalize_total(adata_ct, target_sum=1e4)
sc.pp.log1p(adata_ct)
sc.pp.highly_variable_genes(adata_ct, min_mean=0.0125, max_mean=3, min_disp=0.5)
# adata_ct.raw = adata_ct # unnecessary step
adata_ct = adata_ct[:, adata_ct.var.highly_variable]

## save the PCA covariates without Harmmony
sc.pp.regress_out(adata_ct, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata_ct, max_value=10)
sc.tl.pca(adata_ct, svd_solver="arpack")
df_raw_pcs = pd.DataFrame(adata_ct.obsm["X_pca"])
# remove batch automatically added at the end by integration
df_raw_pcs.index = [re.sub(r"-[0-9]+$", "", cell) for cell in adata_ct.obs.index]
df_raw_pcs.columns = [f"PC{i+1}" for i in df_raw_pcs.columns]

df_raw_pcs.to_csv(pcs_out_file)
###

# sc.pp.regress_out(adata_ct, ["total_counts", "pct_counts_mt"])
# sc.pp.scale(adata_ct, max_value=10)

# # integration
# sc.tl.pca(adata_ct, svd_solver="arpack")
# sce.pp.harmony_integrate(adata_ct, "sequencing_library")
# df_harmony_pcs = pd.DataFrame(adata_ct.obsm["X_pca_harmony"])
# # remove batch automatically added at the end by integration
# df_harmony_pcs.index = [re.sub(r"-[0-9]+$", "", cell) for cell in adata_ct.obs.index]
# df_harmony_pcs.columns = [f"harmony_PC{i+1}" for i in df_harmony_pcs.columns]

# # save PCs
# df_harmony_pcs.to_csv(pcs_out_file)
