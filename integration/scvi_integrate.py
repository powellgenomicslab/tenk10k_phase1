import glob
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

# CellRanger files
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
files = glob.glob(cellranger_dir+"S*")

for file in files:
    sample = file.replace(cellranger_dir, "")
    adata=sc.read_10x_h5(file+"/outs/filtered_feature_bc_matrix.h5")
    adata.var_names_make_unique()
    adata.obs['sample'] = sample
    adata.obs['barcode'] = adata.obs.index
    adata.obs.index = sample + "_" + adata.obs['barcode']
    datasets.append(adata)
    
adata = datasets[0].concatenate(*datasets[1:])
sc.tl.pca(adata, svd_solver='arpack')
scvi.model.SCVI.setup_anndata(adata, layer = "counts",
                              categorical_covariate_keys=["sample"],
                              continuous_covariate_keys=['pct_counts_mt', 'total_counts'])
model = scvi.model.SCVI(adata)
model.train()
adata.obsm['X_scVI'] = model.get_latent_representation()
adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)
sc.pp.neighbors(adata, use_rep = 'X_scVI')
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.5)
sc.pl.umap(adata, color = ['leiden', 'Sample'], frameon = False)

