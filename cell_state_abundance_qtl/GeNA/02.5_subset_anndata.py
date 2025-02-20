import sys
import numpy as np
import pandas as pd
import multianndata as mad
import cna
import scanpy as sc
import matplotlib.pyplot as plt

np.random.seed(0)

# get cell type from cmdline arg 1
celltype = sys.argv[1]
resolution = sys.argv[2]

# subset the object and save

print(f"Subsetting {celltype} from master Anndata...")

# constants
outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

# read in the latest tenk cohort
adata = sc.read(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/298_libraries_concatenated_filtered.h5ad",
    cache=True,
)

major_cell_type_mapping = {
    "B_intermediate": "B",
    "B_memory": "B",
    "B_naive": "B",
    "Plasmablast": "B",
    "NK": "NK",
    "NK_CD56bright": "NK",
    "NK_Proliferating": "NK",
    "CD8_Naive": "CD8_T",
    "CD8_Proliferating": "CD8_T",
    "CD8_TCM": "CD8_T",
    "CD8_TEM": "CD8_T",
    "CD4_CTL": "CD4_T",
    "CD4_Naive": "CD4_T",
    "CD4_Proliferating": "CD4_T",
    "CD4_TCM": "CD4_T",
    "CD4_TEM": "CD4_T",
    "Treg": "CD4_T",
    "dnT": "Unconventional_T",
    "gdT": "Unconventional_T",
    "ILC": "Unconventional_T",
    "MAIT": "Unconventional_T",
    "pDC": "Dendritic",
    "cDC1": "Dendritic",
    "cDC2": "Dendritic",
    "ASDC": "Dendritic",
    "CD14_Mono": "Monocyte",
    "CD16_Mono": "Monocyte",
    "HSPC": "Other",
}

# add major cell type to single cell metadata

adata.obs["major_cell_type"] = [
    major_cell_type_mapping[ct] for ct in adata.obs["wg2_scpred_prediction"]
]

if resolution == "major_cell_types":
    # subset to only target cell type
    adata = adata[adata.obs["major_cell_type"] == celltype]

elif resolution == "minor_cell_types":
    # subset to only target cell type
    adata = adata[adata.obs["wg2_scpred_prediction"] == celltype]

# save counts
adata.layers["counts"] = adata.X.copy()
adata.write(f"{outdir}/data/h5/{resolution}/{celltype}_scanpy.h5ad")
