import glob
import bbknn
import time
import numpy as np
import pandas as pd
import scanpy as sc


# CellRanger files
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
files = glob.glob(cellranger_dir+"S*")

numbers_to_combine = [2,5,10,20,40,64]
data = {"nsamples": numbers_to_combine, "ncells": 0, "ngenes": 0, "bbknn_time": 0,}
summary_df = pd.DataFrame(data)
summary_filename = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/integration/bbknn_times_summary.csv"

for num in numbers_to_combine:
    datasets = []
    for file in files[0:num]:
        sample = file.replace(cellranger_dir, "")
        adata=sc.read_10x_h5(file+"/outs/filtered_feature_bc_matrix.h5")
        adata.var_names_make_unique()
        adata.obs['sample'] = sample
        adata.obs['barcode'] = adata.obs.index
        adata.obs.index = sample + "_" + adata.obs['barcode']
        datasets.append(adata)
    rowIndex = summary_df[summary_df['nsamples']==num].index[0]
    start_time = time.time()
    adata_all = datasets[0].concatenate(*datasets[1:])
    summary_df.at[rowIndex,'concat_time'] = time.time() - start_time
    summary_df.at[rowIndex,'ncells'] = adata_all.shape[0]
    summary_df.at[rowIndex,'ngenes'] = adata_all.shape[1]
    start_time = time.time()
    sc.tl.pca(adata_all, svd_solver='arpack')
    summary_df.at[rowIndex,'pca_time'] = time.time() - start_time
    start_time = time.time()
    bbknn.bbknn(adata_all, batch_key='sample')
    summary_df.at[rowIndex,'bbknn_time'] = time.time() - start_time
    summary_df.to_csv(summary_filename)

print(summary_df.head())

