import sys
import pandas as pd
import scanpy as sc

celltype = sys.argv[1]

# Specify which directory the files generated here will be saved to
output_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/integrated_objects/cpg_anndata/"

# Specify which files this script will generate
# CSV for expression covariates (principal components after Harmony correction)
pcs_out_file = f'{output_dir}/{celltype}_expression_pcs.csv'

# Load integrated AnnData object
input_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/integrated_objects/"
input_file = f"{out_dir}224_libraries/concatenated_gene_info.csv"

adata = sc.read(input_file)

# Extract cell type specific expression
adata_ct = adata[adata.obs['wg2_scpred_prediction'] == celltype]

# normalisation (per cell type)
sc.pp.normalize_total(adata_ct, target_sum=1e4)
sc.pp.log1p(adata_ct)
sc.pp.highly_variable_genes(adata_ct, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_ct.raw = adata_ct
adata_ct = adata_ct[:, adata_ct.var.highly_variable]
sc.pp.regress_out(adata_ct, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_ct, max_value=10)

# integration
sc.tl.pca(adata_ct, svd_solver='arpack')
df_harmony_pcs = pd.DataFrame(adata_ct.obsm['X_pca_harmony'])
df_harmony_pcs.index = adata_ct.obs.index
df_harmony_pcs.columns = [f"harmony_PC{i+1}" for i in df_harmony_pcs.columns]

# save PCs
df_harmony_pcs.to_csv(pcs_out_file)
