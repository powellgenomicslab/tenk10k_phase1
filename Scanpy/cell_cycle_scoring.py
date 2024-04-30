import scanpy as sc
import pandas as pd

adata = sc.read('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info_harmony_umap_notsne.h5ad')

# ribosomal and haemoglobin genes qc metrics 
# adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL")) # ribosomal genes
# adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))# hemoglobin genes.
# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, log1p=False, inplace=True)

# cell cycle scoring
adata = sc.tl.score_genes_cell_cycle(adata)


adata.obs[['total_counts_mt', 'pct_counts_mt', 'total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb',
       'pct_counts_hb']].to_csv('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_cell_cycle.csv')
