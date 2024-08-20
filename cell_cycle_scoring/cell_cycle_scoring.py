import scanpy as sc
import pandas as pd
import re

adata = sc.read('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info_harmony_umap_notsne.h5ad')

# TODO: [done] remove these columns somewhere earlier in the pipeline 
# it seems like they are getting duplicated when the anndata object are concatenated, we get one of each qc column for each pool resulting in +960 columns
# columns_to_drop = [col for col in adata.var.columns if re.match(r'(n_cells_by_counts-\d+|mean_counts-\d+|pct_dropout_by_counts-\d+|total_counts-\d+)', col)]
# adata.var.drop(columns=columns_to_drop, inplace=True)

# TODO: run this per-library, before qc, output to csv and read in the results in the scanpy combiner script. 

# cell cycle gene list from Tirosh et al 2015, downloaded from https://raw.githubusercontent.com/scverse/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
cell_cycle_genes = [x.strip() for x in open('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/regev_lab_cell_cycle_genes.txt')]

# split into s, g2m phase genes, and convert to ENSG id's
# NOTE: only scoring 9 s genes and 6 g2m genes as only these are included in the 3343 highly variable genes...
s_genes = adata.var['gene_ids'][adata.var['gene_name'].isin(cell_cycle_genes[:43])]
g2m_genes = adata.var['gene_ids'][adata.var['gene_name'].isin(cell_cycle_genes[43:])]

# cell cycle scoring
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

adata.obs[['S_score', 'G2M_score', 'phase']].to_csv('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/240_libraries_cellcyle_phase.csv')

# TODO: Might be worth performing the cell cycle scoring using all the cell cycle genes rather than the tiny subset in the HVG's
