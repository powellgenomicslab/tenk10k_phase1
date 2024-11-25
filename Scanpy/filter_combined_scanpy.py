import scanpy as sc
import pandas as pd
import scanpy.external as sce

sc.settings.verbosity = 4  # show logging output

# Output directory
out_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/"

adata = sc.read(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/concatenated_gene_info.h5ad"
)

# remove platelet, erythrocyte, and doublets as predicted by scPred
adata = adata[
    ~adata.obs["wg2_scpred_prediction"].isin(["Platelet", "Eryth", "Doublet"])
]


# Apply qc thresholds to remove cells with low number of counts and low number of genes detected
sc.pp.filter_cells(adata, min_counts=800)
sc.pp.filter_cells(adata, min_genes=600)

# save the filtered object
adata.write(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_concatenated_filtered.h5ad"
)
