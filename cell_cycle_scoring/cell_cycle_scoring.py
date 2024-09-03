import sys
import glob
import re
import scanpy as sc
import pandas as pd

from cellbender.remove_background.downstream import anndata_from_h5

# adata = sc.read('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info_harmony_umap_notsne.h5ad')

# TODO:
# [] run this for all libraries
# [] read in the results in the scanpy combiner script adding them to each object and the combined final scanpy object
# [x] use annotation to convert ensg to gene symbols?
# [x] run using all genes
# [x] output to csv and
# [x] var names not unique warning??

i = int(sys.argv[1])

# input
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"

# output
output_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/per_library_phases"

samples = glob.glob(cellranger_dir + "S*")
sample = samples[i - 1]  # mismatch in index between bash and python
sample = sample.replace(cellranger_dir, "")

# create anndata object for one library
filtered_matrix = f"{cellranger_dir}/{sample}/outs/filtered_feature_bc_matrix.h5"

print(f"reading in filtered matrix {filtered_matrix}...")

adata = sc.read_10x_h5(filtered_matrix)

print("filtering cells...")

# I think I wont do any gene filtering to be consistent with what Anna did with the seurat objects - use all cells and all genes
# cell cycle gene list from Tirosh et al 2015, downloaded from https://raw.githubusercontent.com/scverse/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
cell_cycle_genes = [
    x.strip()
    for x in open(
        "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/regev_lab_cell_cycle_genes.txt"
    )
]

# add gene annotations
gencode_text_file = "/directflow/SCCGGroupShare/projects/anncuo/reference_data/gencode.v44.basic.annotation_df.txt"
gene_info = pd.read_csv(gencode_text_file)

adata.var.index = [gene for gene in adata.var["gene_ids"]]
genes_cellranger = adata.var.index
# Remove gene version
gene_info.index = [gene.split(".")[0] for gene in gene_info["gene_id"]]
# Only consider genes in AnnData object
gene_info_df = gene_info[gene_info.index.isin(genes_cellranger)]
adata.var = pd.concat([adata.var, gene_info_df], axis=1)


# split into s, g2m phase genes, and convert to ENSG id's
# NOTE: only scoring 9 s genes and 6 g2m genes as only these are included in the 3343 highly variable genes...
s_genes = adata.var["gene_ids"][adata.var["gene_name"].isin(cell_cycle_genes[:43])]
g2m_genes = adata.var["gene_ids"][adata.var["gene_name"].isin(cell_cycle_genes[43:])]

# s_genes = cell_cycle_genes[:43]
# g2m_genes = cell_cycle_genes[43:]

# 42 S genes and 52 G2M genes getting used for this sample
print(f"Running cell cycle scoring with:\ns_genes:\n{s_genes}\ng2m_genes:\n{g2m_genes}")

# cell cycle scoring
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
# remove batch from cell index (handling cases of a hyphen in the sample name)
adata.obs.index = [re.sub(r"-[0-9]+$", "", cell) for cell in adata.obs.index]

# write output csv
adata.obs[["S_score", "G2M_score", "phase"]].to_csv(
    f"{output_dir}/{sample}_cell_cycle_scoring.csv"
)
