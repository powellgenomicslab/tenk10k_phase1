# mm activate cellregmap_notebook_bb

import scanpy as sc
import pandas as pd
from matplotlib import rcParams, rc

sc.settings.set_figure_params(dpi=400, dpi_save=500)
sc.settings.figdir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/"

adata = sc.read(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_concatenated_harmony_filtered_min1000genes.h5ad"
)

# set plot order
celltype_order = [
    "CD4 TEM",
    "CD4 CTL",
    "Treg",
    "CD4 TCM",
    "CD4 Naive",
    "CD8 TCM",
    "CD8 Naive",
    "CD8 TEM",
    "dnT",
    "gdT",
    "MAIT",
    "ILC",
    "NK",
    "NK CD56bright",
    "B naive",
    "B intermediate",
    "B memory",
    "Plasmablast",
    "CD14 Mono",
    "CD16 Mono",
    "cDC2",
    "pDC",
    "cDC1",
    "ASDC",
    "HSPC",
    "CD4 Proliferating",
    "CD8 Proliferating",
    "NK Proliferating",
]

# nicer formatting of cell names
adata.obs["Cell type"] = adata.obs["wg2_scpred_prediction"].str.replace("_", " ")
adata.obs["Cell type"] = adata.obs["Cell type"].astype("category")

# color palette tenk10k
palette_df = pd.read_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/colour_palette_table.tsv",
    sep="\t",
)
palette_dict = dict(zip(palette_df["cell_type"], palette_df["color"]))
rc("font", **{"family": "sans-serif", "sans-serif": ["Helvetica"]})
rcParams["figure.figsize"] = (8, 3)
rcParams["font.size"] = 22

sc.pl.violin(
    adata,
    "ENSG00000090382",
    groupby="Cell type",
    palette=palette_dict,
    size=0,
    order=palette_df["cell_type"],
    rotation=90,
    stripplot=False,
    save="lyz_vln.png",
)

sc.pl.violin(
    adata,
    "ENSG00000090382",
    groupby="Cell type",
    palette=palette_dict,
    size=0,
    order=palette_df["cell_type"],
    rotation=90,
    stripplot=False,
    save="lyz_vln.pdf",
)
