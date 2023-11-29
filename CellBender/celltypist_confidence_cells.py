import glob
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
from cellbender.remove_background.downstream import anndata_from_h5

# Download celltypist models
models.download_models(force_update = True)
model = models.Model.load(model = 'Immune_All_Low.pkl')

# Define results directories
cellbender_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/cellbender_output/"
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"

def run_celltypist(adata):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)
    adata_predicted = predictions.to_adata()
    return adata_predicted

samples = [s.replace(cellranger_dir, "") for s in glob.glob(cellranger_dir+"/S*")]

data = {"samples": samples, "cellranger_ncells": 0, "cellranger_avg_celltypist_confidence": 0.0,
       "cellbender_09_ncells": 0, "cellbender_09_avg_celltypist_confidence": 0.0,
        "cellbender_05_ncells": 0, "cellbender_05_avg_celltypist_confidence": 0.0,}
summary_df = pd.DataFrame(data)
summary_filename = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/cellbender/cells_set_celltypist_confidence_summary.csv"

for sample in samples:
    print(sample)
    # define files
    cellbender_file = cellbender_dir + sample + "/cellbender_output.h5"
    cellranger_raw_file = cellranger_dir + sample + "/outs/raw_feature_bc_matrix.h5"
    cellranger_filt_file = cellranger_dir + sample + "/outs/filtered_feature_bc_matrix.h5"
    # all cells from filtered object
    adata_filt = sc.read_10x_h5(cellranger_filt_file)
    cellranger_cells = adata_filt.obs.index
    # load cellbender output
    cellbender_adata = anndata_from_h5(cellbender_file)
    cellbender_cells_09 = cellbender_adata[cellbender_adata.obs['cell_probability'] > 0.9].obs.index
    cellbender_cells_05 = cellbender_adata[cellbender_adata.obs['cell_probability'] > 0.5].obs.index
    # subset objects to different cells
    adata = sc.read_10x_h5(cellranger_raw_file)
    adata_cellranger = adata[adata.obs.index.isin(cellranger_cells)]
    adata_cellbender_09 = adata[adata.obs.index.isin(cellbender_cells_09)]
    adata_cellbender_05 = adata[adata.obs.index.isin(cellbender_cells_05)]
    adata_cellranger_predicted = run_celltypist(adata_cellranger)
    adata_cellbender_09_predicted = run_celltypist(adata_cellbender_09)
    adata_cellbender_05_predicted = run_celltypist(adata_cellbender_05)
    # enter values into summary table
    rowIndex = summary_df[summary_df['samples']==sample].index[0]
    summary_df.at[rowIndex,'cellranger_ncells'] = adata_cellranger_predicted.obs.shape[0]
    summary_df.at[rowIndex,'cellbender_09_ncells'] = adata_cellbender_09_predicted.obs.shape[0]
    summary_df.at[rowIndex,'cellbender_05_ncells'] = adata_cellbender_05_predicted.obs.shape[0]
    summary_df.at[rowIndex,'cellranger_avg_celltypist_confidence'] = adata_cellranger_predicted.obs['conf_score'].mean()
    summary_df.at[rowIndex,'cellbender_09_avg_celltypist_confidence'] = adata_cellbender_09_predicted.obs['conf_score'].mean()
    summary_df.at[rowIndex,'cellbender_05_avg_celltypist_confidence'] = adata_cellbender_05_predicted.obs['conf_score'].mean()
    summary_df.to_csv(summary_filename)
