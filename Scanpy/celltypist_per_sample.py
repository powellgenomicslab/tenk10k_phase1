import glob
import scanpy as sc

import celltypist
from celltypist import models

# CellRanger files
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
files = glob.glob(cellranger_dir+"S*")

myfile = files[0]

# load adata
adata=sc.read_10x_h5(myfile+"/outs/filtered_feature_bc_matrix.h5")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# confirm normalisation of AnnData
adata.X.expm1().sum(axis = 1)

# download all models 
models.download_models(force_update = True)

# specify immune (sub)cell type model
model = models.Model.load(model = 'Immune_All_Low.pkl')

# predict cell types and save as adata
predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)
adata_predicted = predictions.to_adata()

# plot UMAP with cell types annotations
sc.pl.umap(adata_predicted, color = 'majority_voting')
