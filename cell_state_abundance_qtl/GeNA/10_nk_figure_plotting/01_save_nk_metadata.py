import scanpy as sc
import pandas as pd

# export metadata from lawrence's NK cell analysis

nk_path = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/NK_cellstate/ALL_TenK10K_NK_rerun_maxgs400_ep50_beta1_annealKLTrue_lr0.02_batchsize32_umap.h5ad"

adata = sc.read(nk_path, cache=True)
adata.obs.columns

umap_df = pd.DataFrame(
    adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
).join(adata.obs.loc[:, ["celltype"]])

umap_df.to_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/NK_cellstate/NK_cellstate_umap_subtypes.csv"
)
