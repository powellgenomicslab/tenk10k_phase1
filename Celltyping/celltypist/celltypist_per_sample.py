import glob
import os
import sys

import scanpy as sc

import celltypist
from celltypist import models

i = int(sys.argv[1])
SEQ_DATE=str(sys.argv[2])

# CellRanger files 

# 64 samples from 231013
# cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"

# 24 samples from 231213
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231213_tenk10k_gencode44/cellranger_outs/"

# 41 samples from 231214
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231214_tenk10k_gencode44/cellranger_outs/"

# 18 samples from 240108
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240108_tenk10k_gencode44/cellranger_outs/"

# 18 samples from 240112
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240112_tenk10k_gencode44/cellranger_outs/"

# 25 samples from 240115
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240115_tenk10k_gencode44/cellranger_outs/"

# 17 samples from 240116
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240116_tenk10k_gencode44/cellranger_outs/"

# 17 samples from 240119
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240119_tenk10k_gencode44/cellranger_outs/"

# 16 samples from 240214
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240214_tenk10k_gencode44/cellranger_outs/"

# 17 samples from 240223

# use sequencing date supplied via arg parser rather than manually updating this script 
cellranger_dir = f"/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/{SEQ_DATE}_tenk10k_gencode44/cellranger_outs/"

# Output directory
output_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/celltypist/output/"

samples = glob.glob(cellranger_dir+"S*")
# mismatch in index between bash and python
sample = samples[i-1]
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

# umap
sc.tl.umap(adata_predicted)

# save
adata_predicted.write(predicted_filename)
