import glob
import re
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce

# Output directory
out_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/"

# Filtered object directory
scanpy_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/scanpy_objects_w_metadata/"

# extract all files
scanpy_files = glob.glob(scanpy_dir + "S*")

# combine anndata objects by samples
datasets = []
for file in scanpy_files:
    adata = sc.read(file)

    # TODO: [DONE] REMOVE QC columns before merging, otherwise the are duplicated for each pool and make the final object way bigger
    # columns_to_drop = [col for col in adata.var.columns if re.match(r'(n_cells_by_counts-\d+|mean_counts-\d+|pct_dropout_by_counts-\d+|total_counts-\d+)', col)]
    # adata.var.drop(columns=columns_to_drop, inplace=True)

    datasets.append(adata)
adata = datasets[0].concatenate(*datasets[1:])

# Load gencode (v44) gene info file
gencode_text_file = (
    # "/share/ScratchGeneral/anncuo/reference_data/gencode.v44.basic.annotation_df.txt"
    "/directflow/SCCGGroupShare/projects/anncuo/reference_data/gencode.v44.basic.annotation_df.txt"
)
gene_info = pd.read_csv(gencode_text_file)

# Add gene info to AnnData
# set gene ids as indices
adata.var.index = [gene for gene in adata.var["gene_ids"]]
genes_cellranger = adata.var.index
# Remove gene version
gene_info.index = [gene.split(".")[0] for gene in gene_info["gene_id"]]
# Only consider genes in AnnData object
gene_info_df = gene_info[gene_info.index.isin(genes_cellranger)]
# Add info
adata.var = pd.concat([adata.var, gene_info_df], axis=1)

# remove batch from cell index (handling cases of a hyphen in the sample name)
adata.obs.index = [re.sub(r"-[0-9]+$", "", cell) for cell in adata.obs.index]

# write
out_file = f"{out_dir}300_libraries/concatenated_gene_info.h5ad"
adata.write(out_file)
