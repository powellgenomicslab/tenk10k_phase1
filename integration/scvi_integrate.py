import glob
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix

# Filtered object directory
scanpy_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/filtered_scanpy_objects/"

# load adata
scanpy_files = glob.glob(scanpy_dir+"S*")

datasets = []
for file in scanpy_files:
    adata = sc.read(file)
    datasets.append(adata)
    
adata = datasets[0].concatenate(*datasets[1:])

adata.X = csr_matrix(adata.X)

sc.pp.filter_genes(adata, min_cells = 100)

adata.layers['counts'] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)
adata.raw = adata

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
sc.pl.umap(adata, color = ['leiden', 'sample'], frameon = False)
adata.uns['scvi_markers'] = markers_scvi
adata.uns['markers'] = markers
adata.write_h5ad('integrated.h5ad')
model.save('model.model')

