import glob
import re
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from cellbender.remove_background.downstream import anndata_from_h5

# CellRanger files
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
cellranger_files = glob.glob(cellranger_dir+"S*")

# Cellbender files
cellbender_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/cellbender_output_smaller_learning_rate/"

# Demuxafy files (combined results)
demuxafy_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/combined_output_scds_scdblfinder_vireo/"

# Output directory
output_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/filtered_scanpy_objects/"

for sample_path in cellranger_files:
    sample = sample_path.replace(cellranger_dir,"")
    # Load Cellranger counts
    adata=sc.read_10x_h5(sample_path+"/outs/filtered_feature_bc_matrix.h5")
    # Load Cellbender file
    cellbender_file = cellbender_dir + sample + "/" + sample + "cellbender_output.h5"
    cellbender_adata = anndata_from_h5(cellbender_file)
    # subset to cells estimated by Cellranger
    cells_cellranger = adata.obs.index
    cellbender_df = cellbender_adata.obs[cellbender_adata.obs.index.isin(cells_cellranger)]
    # rename columns to specify where the info comes from
    cellbender_df.columns = ["cellbender_" + i for i in cellbender_df.columns]
    # add cellbender info to adata obs
    adata.obs = pd.concat([adata.obs,cellbender_df], axis=1)
    # get combined demultiplexing + doublet info
    demuxafy_file = demuxafy_dir + sample + "/combined_results_w_combined_assignments.tsv"
    demuxafy_df = pd.read_csv(demuxafy_file, sep="\t")
    demuxafy_df.index = demuxafy_df['Barcode']
    demuxafy_df = demuxafy_df.drop(columns=['Barcode'])
    # add this info to the obs as well
    adata.obs = pd.concat([adata.obs,demuxafy_df], axis=1)
    # filter data to cells that are singlets and assigned to an individual 
    adata = adata[~adata.obs['MajoritySinglet_Individual_Assignment'].isin(["unassigned","doublet"])]
    # preprocessing
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
    adata.obs["individual"] = sample
    donors = adata.obs['MajoritySinglet_Individual_Assignment']
    # add samples ID to donors that are just added by vireo to make them unique
    donor_regex = re.compile(r'donor[0-9]+')
    adata.obs["individual"] = [f"{i}_{sample}" if re.match(donor_regex, i) else i for i in donors]
    adata.write(output_dir+sample+"_w_cellbender_no_doublet_no_unassigned_basic_pp.h5ad")
