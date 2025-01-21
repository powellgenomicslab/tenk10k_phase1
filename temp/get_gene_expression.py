# corresponding metadata: "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_subset_filtered.csv"

import scanpy as sc

filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_concatenated_harmony_filtered.h5ad"
adata = sc.read(filename)
adata = adata.raw.to_adata()  # get the raw counts

df = sc.get.obs_df(adata, keys="ENSG00000089041")
df.to_csv("/directflow/SCCGGroupShare/projects/blabow/temp/P2RX7_counts.csv")
