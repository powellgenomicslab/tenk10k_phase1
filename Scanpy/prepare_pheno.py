import os
import sys
import pandas as pd
import scanpy as sc
import re

celltype = sys.argv[1]
chromosome = sys.argv[2]

print(celltype)
print(chromosome)

# integrated_objects_dir='/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries'
integrated_objects_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries"
# Specify which files this script will generate
output_dir = f"{integrated_objects_dir}/cpg_anndata_filtered"

if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)

# Specify which directory the files generated here will be saved to
# AnnData object: expression + gene info + batch info (batch=sequencing library)
# ct = celltype.replace(' ','_') # remove spaces from cell type names
adata_out_file = f"{output_dir}/{celltype}_chr{chromosome}.h5ad"
# if os.path.exists(adata_out_file):
#   sys.exit('File already exists!')

# Load integrated AnnData object
# input_file = f'{integrated_objects_dir}/240_libraries_concatenated_harmony_leiden_filtered_reanalysed.h5ad'
input_file = f"{integrated_objects_dir}/300_libraries_concatenated_filtered.h5ad"
adata = sc.read(input_file)

# Extract cell type and chromosome specific expression
adata_ct = adata[adata.obs["wg2_scpred_prediction"] == celltype]
adata_ct_chr = adata_ct[:, adata_ct.var["chr"] == f"chr{chromosome}"]
adata_ct_chr.index = [
    re.sub(r"-[0-9]+$", "", cell) for cell in adata_ct_chr.obs.index
]  # now this happens in combine_files_add_gene_info.py can maybe remove this line

# adata_ct_chr.raw = adata_ct_chr.copy()

# write
adata_ct_chr.write(adata_out_file)
