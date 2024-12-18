import scanpy as sc
import pandas as pd
import scanpy.external as sce

sc.settings.verbosity = 4  # show logging output

# Output directory
out_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/"

adata = sc.read(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_concatenated_filtered.h5ad"
)

# get a subset of the metadata for UMAP QC plots
adata.obs.to_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_filtered_min1000genes.csv"
)
