import os
import sys
import pandas as pd
import scanpy as sc

celltype = sys.argv[1]
chromosome = sys.argv[2]

print(celltype)
print(chromosome)

integrated_objects_dir='/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects'
# Specify which files this script will generate
output_dir = f'{integrated_objects_dir}/cpg_anndata'
# Specify which directory the files generated here will be saved to
# AnnData object: expression + gene info + batch info (batch=sequencing library)
# ct = celltype.replace(' ','_') # remove spaces from cell type names
adata_out_file = f'{output_dir}/{celltype}_chr{chromosome}.h5ad'
if os.path.exists(adata_out_file):
  sys.exit('File already exists!')

# Load integrated AnnData object
input_file = f'{integrated_objects_dir}/240_libraries_concatenated_gene_info.h5ad'
adata = sc.read(input_file)

# Extract cell type and chromosome specific expression
adata_ct = adata[adata.obs['wg2_scpred_prediction'] == celltype]
adata_ct_chr = adata_ct[:, adata_ct.var['chr'] == f'chr{chromosome}']
adata_ct_chr.obs['cell'] = adata_ct_chr.obs.index

# write
adata_ct_chr.write(adata_out_file)
