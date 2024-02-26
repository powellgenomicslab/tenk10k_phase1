import sys
import pandas as pd
import scanpy as sc

celltype = sys.argv[1]
chromosome = sys.argv[2]

print(celltype)
print(chromosome)

# Specify which files this script will generate
output_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/integrated_objects/cpg_anndata"
# Specify which directory the files generated here will be saved to
# AnnData object: expression + gene info + batch info (batch=sequencing library)
ct = celltype.replace(" ","_") # remove spaces from cell type names
adata_out_file = f'{output_dir}/{ct}_chr{chromosome}.h5ad'

print(adata_out_file)

# Load integrated AnnData object
input_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/integrated_objects/"
input_file = f"{input_dir}224_libraries/concatenated_gene_info_donor_info.h5ad"
adata = sc.read(input_file)

# Extract cell type and chromosome specific expression
adata_ct = adata[adata.obs['wg2_scpred_prediction'] == celltype]
adata_ct_chr = adata_ct[:, adata_ct.var["chr"] == f'chr{chromosome}']
adata_ct_chr.obs['cell'] = adata_ct_chr.obs.index

# write
adata_ct_chr.write(adata_out_file)
