import glob
import os
import re
import sys

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from cellbender.remove_background.downstream import anndata_from_h5

# TODO:
# [] add in cell cycle phases csv info to scanpy objects (/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/per_library_phases/)

# index provided, running one sequencing library at a time
i = int(sys.argv[1])
# seq_date = str(sys.argv[2]) # sequencing date is 2nd argument ?

# CellRanger files

# 64 samples from 231013
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231013_tenk10k_gencode44/cellranger_outs/"
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


# set to new for 240214 onwards; note these may get moved in the future
source_data_location = "new"  # set to 'old' or 'new' to read from blake / annna's directory structure  # submit number of jobs = number of files in the cellranger dir

# source_data_location = "old"

# if running on 'new' data also specify the seq_date

# seq_date = "240223"
seq_date = "240501"
# record of previous runs:
# seq_date = "240214"
# seq_date = "240524"

# get cellranger outputs
cellranger_dir = f"/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/{seq_date}_tenk10k_gencode44/cellranger_outs/"
# older cellranger outs moved to here:
# cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"

print(cellranger_dir)


# print(seq_date)

if source_data_location == "new":

    # get sample names from the cellranger_outs txt file
    sample_names_file = open(
        f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cellranger_outs_{seq_date}.txt",
        "r",
    )
    sample_names_list = sample_names_file.read().split("\n")
    samples = [
        s for s in sample_names_list if s
    ]  # remove the blank value at the end of the list
    sample = samples[i - 1]  # mismatch in index between bash and python

    # Cellbender files
    cellbender_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cellbender/output/smaller_learning_rate/"
    cellbender_file = cellbender_dir + sample + "/cellbender_output.h5"

    # Demuxafy files (combined results, vireo w/o cb)
    demuxafy_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/combine_results/output/combined_output_scds_scdblfinder_vireo_no_cb/"

    # Celltypist files
    celltypist_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/celltypist/output/"

    # scPred files
    scpred_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scpred/"


elif source_data_location == "old":

    samples = glob.glob(cellranger_dir + "S*")
    # # mismatch in index between bash and python
    sample = samples[i - 1]
    sample = sample.replace(cellranger_dir, "")

    # Cellbender files
    cellbender_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/cellbender_output_smaller_learning_rate/"

    # some of these have the sample names appended to the start
    cellbender_file = cellbender_dir + sample + "/" + sample + "cellbender_output.h5"

    if not os.path.exists(cellbender_file):
        cellbender_file = cellbender_dir + sample + "/" + "cellbender_output.h5"

    # Demuxafy files (combined results, vireo w/o cb)
    demuxafy_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/combined_output_scds_scdblfinder_vireo_no_cb/"

    # Celltypist files
    celltypist_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/celltypist/"

    # scPred files
    scpred_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/"


# Output directory
# old outdir used for saige-qtl pipeline (with WGS qc fails, odd celltype composition donors removed)
# output_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/scanpy_objects_w_metadata/"

# output dir with all donors except those that withdrew consent
output_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/scanpy_objects_w_metadata/all_donors/"

print(f"Combining metadata for {sample}")

output_filename = output_dir + sample + "_w_metadata_donor_info.h5ad"
# if os.path.exists(output_filename):
#   sys.exit("File already exists!")

# Load Cellranger counts
# filtered_matrix = cellranger_dir+sample+"/outs/filtered_feature_bc_matrix.h5"
# filtered_matrix = cellranger_dir+sample+"/cellranger_count/"+sample+"/outs/filtered_feature_bc_matrix.h5"

filtered_matrix = f"{cellranger_dir}/{sample}/outs/filtered_feature_bc_matrix.h5"

if not os.path.exists(filtered_matrix):
    filtered_matrix = f"{cellranger_dir}/{sample}/cellranger_count/{sample}/outs/filtered_feature_bc_matrix.h5"  # use this location if data not moved yet

adata = sc.read_10x_h5(filtered_matrix)

# Load Cellbender file
cellbender_adata = anndata_from_h5(cellbender_file)

# subset to cells estimated by Cellranger
cells_cellranger = adata.obs.index
cellbender_df = cellbender_adata.obs[cellbender_adata.obs.index.isin(cells_cellranger)]

# rename columns to specify where the info comes from
cellbender_df.columns = ["cellbender_" + i for i in cellbender_df.columns]

# add cellbender info to adata obs
adata.obs = pd.concat([adata.obs, cellbender_df], axis=1)

# add cell typing info
# celltypist
celltypist_file = celltypist_dir + sample + "_celltypist_predicted.h5"
celltypist_adata = sc.read(celltypist_file)
celltypist_df = celltypist_adata.obs[celltypist_adata.obs.index.isin(cells_cellranger)]

# rename columns to specify where the info comes from
celltypist_df.columns = ["celltypist_" + i for i in celltypist_df.columns]

# add celltypist info to adata obs
adata.obs = pd.concat([adata.obs, celltypist_df], axis=1)

# azimuth + hierarchical scPred
scpred_file = scpred_dir + sample + "/combined_metadata.csv"
scpred_df = pd.read_csv(scpred_file, index_col=0)
scpred_df = scpred_df[scpred_df.index.isin(cells_cellranger)]

# rename columns to specify where the info comes from
scpred_df.drop(
    ["orig.ident.1", "nCount_RNA.1", "nFeature_RNA.1", "percent.mt.1"],
    axis=1,
    inplace=True,
)
scpred_df.rename(
    columns={
        "orig.ident": "sample",
        "percent.mt": "percent_mt",
        "predicted.celltype.l2": "azimuth_predicted_celltype_l2",
        "predicted.celltype.l2.score": "azimuth_predicted_celltype_l2_score",
    },
    inplace=True,
)
celltype_columns = ["azimuth_predicted_celltype_l2", "scpred_prediction"]
scpred_df[celltype_columns] = scpred_df[celltype_columns].apply(
    lambda x: x.str.replace(" ", "_")
)
scpred_df.columns = ["wg2_" + i for i in scpred_df.columns]

# add scpred info to adata obs
adata.obs = pd.concat([adata.obs, scpred_df], axis=1)

# get combined demultiplexing + doublet info
demuxafy_file = demuxafy_dir + sample + "/combined_results_w_combined_assignments.tsv"
demuxafy_df = pd.read_csv(demuxafy_file, sep="\t")
demuxafy_df.index = demuxafy_df["Barcode"]
demuxafy_df.drop(columns=["Barcode"], inplace=True)
demuxafy_df = demuxafy_df[demuxafy_df.index.isin(cells_cellranger)]

# add this info to the obs as well
adata.obs = pd.concat([adata.obs, demuxafy_df], axis=1)

# filter data to cells that are singlets and assigned to an individual
adata = adata[
    ~adata.obs["MajoritySinglet_Individual_Assignment"].isin(["unassigned", "doublet"])
]
adata = adata[adata.obs["MajoritySinglet_Individual_Assignment"].notna()]

# preprocessing & QC plotting
adata.var_names_make_unique()
# sc.pp.filter_cells(adata, min_genes=200)
# calculate mitochondrial, ribosomal and hemoglobin gene expression
adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)
plt.savefig(output_dir + "figures/" + sample + "_qc_violins.pdf")

sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
plt.savefig(output_dir + "figures/" + sample + "_tot_counts_pct_mt.pdf")

sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
plt.savefig(output_dir + "figures/" + sample + "_tot_counts_tot_genes.pdf")

# add cell, individual and sequencing library id pre-merging
adata.obs["original_barcode"] = adata.obs.index
adata.obs["new_cell_name"] = [
    i.replace("-1", f"_{sample}") for i in adata.obs["original_barcode"]
]
adata.obs["sequencing_library"] = sample.replace("-", "_")
adata.obs["individual"] = adata.obs["MajoritySinglet_Individual_Assignment"]

donor_regex = re.compile(r"donor[0-9]+")
# remove donors that were not able to be assigned by vireo due to missing genotype data
adata = adata[~adata.obs["individual"].str.contains(donor_regex, regex=True, na=False)]

# add cohort info
seqlib_cohort_map_file = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/sequencing_library_to_cohort_map_300_libs.csv"
seqlib_cohort_map_df = pd.read_csv(seqlib_cohort_map_file)
cohort = seqlib_cohort_map_df[seqlib_cohort_map_df["sequencing_library"] == sample][
    "cohort"
].values[0]
adata.obs["cohort"] = cohort

# add cell cycle scoring
cell_cycle_scoring_csv = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/per_library_phases/output/{sample}_cell_cycle_scoring.csv"
cell_cycle_df = pd.read_csv(cell_cycle_scoring_csv, index_col=0)
adata.obs = adata.obs.join(cell_cycle_df)

# add individual info based on cohort
if cohort == "TOB":
    # map onek1k ids to cpg ids
    df_samples_file = "/directflow/SCCGGroupShare/projects/anncuo/OneK1K_from_ScratchGeneral/scrna-seq_grch38_association_files_OneK1K_CPG_IDs.tsv"
    df_samples = pd.read_csv(df_samples_file, sep="\t")
    df_samples.columns = ["onek1k_id", "cpg_id_old", "tob_id"]
    # the CPG in that file are actually old, so need an extra step to update to the more recent ones (using the same file as below)
    cpg_map_file = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/str_sample-sex-mapping_sample_karyotype_sex_mapping.csv"
    cpg_map_df = pd.read_csv(cpg_map_file)
    cpg_map_df.columns = ["cpg_id", "tob_id", "sex_karyotype"]
    cpg_map_df["tob_id"] = cpg_map_df["tob_id"].str.replace(
        "-PBMC", ""
    )  # remove the trailing string present on only some TOB ID's
    cpg_map_df.drop(columns=["sex_karyotype"], inplace=True)  # drop this column
    df_samples2 = df_samples.merge(cpg_map_df, on="tob_id", how="left")
    df_samples2["individual"] = [
        donor.split("_")[-1] for donor in df_samples2["onek1k_id"]
    ]
    adata.obs = adata.obs.merge(df_samples2, on="individual", how="left")

    adata.obs["onek1k_donor"] = adata.obs["individual"]
    adata.obs["individual"] = adata.obs["cpg_id"]

if cohort == "BioHEART":

    cpg_map_file = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/metadata/bioheart_ct_ids_to_cpg_ids_20241113_n1419.csv"
    cpg_map_df = pd.read_csv(cpg_map_file, index_col=0)
    cpg_map_df.columns = ["ct_id", "individual"]

    adata.obs = adata.obs.merge(cpg_map_df, on="individual", how="left")
    adata.obs["cpg_id"] = adata.obs["individual"]


adata.obs.index = [donor for donor in adata.obs["new_cell_name"]]

# cell type and individual info are key, so drop if NA
adata = adata[adata.obs["individual"].notna()]
adata = adata[adata.obs["wg2_scpred_prediction"].notna()]

# remove data for participants who withdrew consent: "CT_557", "CT_1545", "CT_888"
if cohort == "BioHEART":
    adata = adata[~adata.obs["ct_id"].isin(["CT_557", "CT_1545", "CT_888"])]


# remove the gene-level QC columns as these get duplicated when concatenating the scanpy objects
columns_to_drop = [
    col
    for col in adata.var.columns
    if re.match(
        r"(n_cells_by_counts|mean_counts|pct_dropout_by_counts|total_counts)", col
    )
]
adata.var.drop(columns=columns_to_drop, inplace=True)

# save
adata.write(output_filename)
