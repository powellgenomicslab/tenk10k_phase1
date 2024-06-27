



import scanpy as sc
import pandas as pd
import scanpy.external as sce

sc.settings.verbosity = 4  # show logging output

adata = sc.read('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony.h5ad')

# NOTE: probably will be better to upgrade to newer version of scanpy to use the faster igraph implementation
sc.tl.leiden(
    adata,
)

# NOTE: we previously replaced the PCA with the Batch-corrected Harmony PCA 

adata.write('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony_leiden.h5ad')