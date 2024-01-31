import glob
import os
import re
import sys

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from cellbender.remove_background.downstream import anndata_from_h5

i = int(sys.argv[1])

# CellRanger files 

# 64 samples from 231013
# cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"

# 24 samples from 231213
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231213_tenk10k_gencode44/cellranger_outs/"

# 41 samples from 231214
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231214_tenk10k_gencode44/cellranger_outs/"

# 18 samples from 240108
cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240108_tenk10k_gencode44/cellranger_outs/"

# 18 samples from 240112
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240112_tenk10k_gencode44/cellranger_outs/"

# 25 samples from 240115
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240115_tenk10k_gencode44/cellranger_outs/"

# 17 samples from 240116
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240116_tenk10k_gencode44/cellranger_outs/"

# 17 samples from 240119
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240119_tenk10k_gencode44/cellranger_outs/"

cellranger_files = glob.glob(cellranger_dir+"S*")

samples = glob.glob(cellranger_dir+"S*")
# mismatch in index between bash and python
sample = samples[i-1]
sample = sample.replace(cellranger_dir,"")
print(sample)

# Cellbender files
cellbender_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/cellbender_output_smaller_learning_rate/"

# Demuxafy files (combined results)
demuxafy_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/combined_output_scds_scdblfinder_vireo/"

# Celltypist files
celltypist_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/celltypist/"

# scPred files
scpred_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/"

# Output directory
output_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/filtered_scanpy_objects/"

output_filename = output_dir+sample+"_w_metadata_basic_pp.h5ad"
if os.path.exists(output_filename):
  sys.exit("File already exists!")

# Load Cellranger counts
# filtered_matrix = cellranger_dir+sample+"/outs/filtered_feature_bc_matrix.h5"
filtered_matrix = cellranger_dir+sample+"/cellranger_count/"+sample+"/outs/filtered_feature_bc_matrix.h5"
adata=sc.read_10x_h5(filtered_matrix)

# Load Cellbender file
cellbender_file = cellbender_dir + sample + "/" + sample + "cellbender_output.h5"
# cellbender_file = cellbender_dir + sample + "/cellbender_output.h5"
cellbender_adata = anndata_from_h5(cellbender_file)

# subset to cells estimated by Cellranger
cells_cellranger = adata.obs.index
cellbender_df = cellbender_adata.obs[cellbender_adata.obs.index.isin(cells_cellranger)]

# rename columns to specify where the info comes from
cellbender_df.columns = ["cellbender_" + i for i in cellbender_df.columns]

# add cellbender info to adata obs
adata.obs = pd.concat([adata.obs,cellbender_df], axis=1)

# add cell typing info
# celltypist
celltypist_file = celltypist_dir + sample + "_celltypist_predicted.h5"
celltypist_adata = sc.read(celltypist_file)
celltypist_df = celltypist_adata.obs[celltypist_adata.obs.index.isin(cells_cellranger)]

# rename columns to specify where the info comes from
celltypist_df.columns = ["celltypist_" + i for i in celltypist_df.columns]

# add celltypist info to adata obs
adata.obs = pd.concat([adata.obs,celltypist_df], axis=1)

# azimuth + hierarchical scPred
scpred_file = scpred_dir + sample + "/combined_metadata.csv"
scpred_df = pd.read_csv(scpred_file, index_col=0)
scpred_df = scpred_df[scpred_df.index.isin(cells_cellranger)]

# rename columns to specify where the info comes from
scpred_df.drop(["orig.ident.1","nCount_RNA.1","nFeature_RNA.1","percent.mt.1"], axis=1, inplace=True)
scpred_df.rename(columns={'orig.ident': 'sample', 'percent.mt': 'percent_mt',
                         'predicted.celltype.l2': 'azimuth_predicted_celltype_l2',
                         'predicted.celltype.l2.score': 'azimuth_predicted_celltype_l2_score'}, inplace=True)
scpred_df.columns = ["wg2_" + i for i in scpred_df.columns]

# add scpred info to adata obs
adata.obs = pd.concat([adata.obs,scpred_df], axis=1)

# get combined demultiplexing + doublet info
demuxafy_file = demuxafy_dir + sample + "/combined_results_w_combined_assignments.tsv"
demuxafy_df = pd.read_csv(demuxafy_file, sep="\t")
demuxafy_df.index = demuxafy_df['Barcode']
demuxafy_df.drop(columns=['Barcode'], inplace=True)
demuxafy_df = demuxafy_df[demuxafy_df.index.isin(cells_cellranger)]

# add this info to the obs as well
adata.obs = pd.concat([adata.obs,demuxafy_df], axis=1)

# filter data to cells that are singlets and assigned to an individual 
adata = adata[~adata.obs['MajoritySinglet_Individual_Assignment'].isin(["unassigned","doublet"])]

# preprocessing & QC plotting
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
plt.savefig(output_dir+"figures/"+sample+"_qc_violins.pdf")

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
plt.savefig(output_dir+"figures/"+sample+"_tot_counts_pct_mt.pdf")

sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
plt.savefig(output_dir+"figures/"+sample+"_tot_counts_tot_genes.pdf")

adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
sc.pp.normalize_total(adata, target_sum=1e4) #normalize every cell to 10,000 UMI
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
plt.savefig(output_dir+"figures/"+sample+"_hvgs.pdf")

# add individual and sample id pre-merging
adata.obs["sample"] = sample
adata.obs["individual"] = adata.obs['MajoritySinglet_Individual_Assignment']

# add samples ID to donors that are just added by vireo to make them unique
donor_regex = re.compile(r'donor[0-9]+')
adata.obs["individual"] = [f"{i}_{sample}" if re.match(donor_regex, str(i)) else i for i in adata.obs["individual"]]

# save
adata.write(output_filename)
