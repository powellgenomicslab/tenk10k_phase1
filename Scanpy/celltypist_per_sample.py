import glob
import os
import sys

import scanpy as sc

import celltypist
from celltypist import models

i = int(sys.argv[1])

# CellRanger files 

# 64 samples from 231013
# cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"

# 24 samples from 231213
cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231213_tenk10k_gencode44/cellranger_outs/"

# 41 samples from 231214
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231214_tenk10k_gencode44/cellranger_outs/"


# Ouput directory
output_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/celltypist/"

samples = glob.glob(cellranger_dir+"S*")
# mismatch in index between bash and python
sample = samples[i-1]
print(sample)
sample = sample.replace(cellranger_dir,"")
print(sample)

predicted_filename = output_dir+sample+"_celltypist_predicted.h5"
if os.path.exists(predicted_filename):
  sys.exit("File already exists!")

# load adata
# filtered_matrix = cellranger_dir+sample+"/outs/filtered_feature_bc_matrix.h5"
filtered_matrix = cellranger_dir+sample+"/cellranger_count/"+sample+"/outs/filtered_feature_bc_matrix.h5"
adata=sc.read_10x_h5(filtered_matrix)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# confirm normalisation of AnnData
adata.X.expm1().sum(axis = 1)

# specify immune (sub)cell type model
model = models.Model.load(model = 'Immune_All_Low.pkl')

# predict cell types and save as adata
predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)
adata_predicted = predictions.to_adata()

# plot UMAP with cell types annotations
sc.pl.umap(adata_predicted, color = 'majority_voting')

# save
adata_predicted.write(predicted_filename)
