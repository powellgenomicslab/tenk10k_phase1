import sys
import pandas as pd
import scanpy as sc

celltype = sys.argv[1]
chromosome = sys.argv[2]

# Specify which directory the files generated here will be saved to
input_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/integrated_objects/"
input_file = f"{out_dir}224_libraries/concatenated_gene_info.csv"

# Specify which files this script will generate
output_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/integrated_objects/cpg_anndata/"

# AnnData object: expression + gene info + batch info (batch=sequencing library)
adata_out_file = f'{output_dir}/{celltype}_chr{chromosome}.h5ad'

# CSV for expression covariates (principal components after Harmony correction)
pcs_out_file = f'{out_dir}/{celltype}_expression_pcs.csv'

##############################################################################################

# Load integrated AnnData object
adata = sc.read(input_file)

# Extract cell type and chromosome specific expression
adata_ct = adata[adata.obs['wg2_scpred_prediction'] == celltype]
adata_ct_chr = adata_ct[:, adata_ct.var["chr"] == f'chr{chromosome}']
adata_ct_chr.obs['cell'] = adata_ct_chr.obs.index
adata_ct_chr.write(adata_out_file)


# normalisation
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

# integration
sc.tl.pca(adata_CD4_Naive, svd_solver='arpack')
df_harmony_pcs = pd.DataFrame(adata_CD4_Naive.obsm['X_pca_harmony'])
df_harmony_pcs.index = adata_CD4_Naive.obs.index
df_harmony_pcs.columns = [f"harmony_PC{i+1}" for i in df_harmony_pcs.columns]

# save PCs
df_harmony_pcs.to_csv(pcs_out_file)
