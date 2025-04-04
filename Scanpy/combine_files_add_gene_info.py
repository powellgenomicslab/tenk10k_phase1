import glob
import re
import pandas as pd
import scanpy as sc
import anndata as ad
import scanpy.external as sce

# Output directory
out_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/"

# Filtered object directory
# scanpy_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/scanpy_objects_w_metadata/"

# object directory with all donors
scanpy_dir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/scanpy_objects_w_metadata/all_donors/"

# extract all files
scanpy_files = glob.glob(scanpy_dir + "S*")

# combine anndata objects by samples
datasets = []
for file in scanpy_files:
    adata = sc.read(file)

    # TODO: [DONE] REMOVE QC columns before merging, otherwise the are duplicated for each pool and make the final object way bigger
    # columns_to_drop = [col for col in adata.var.columns if re.match(r'(n_cells_by_counts-\d+|mean_counts-\d+|pct_dropout_by_counts-\d+|total_counts-\d+)', col)]
    # adata.var.drop(columns=columns_to_drop, inplace=True)

    datasets.append(adata)
adata = datasets[0].concatenate(*datasets[1:])

# Load gencode (v44) gene info file
gencode_text_file = (
    # "/share/ScratchGeneral/anncuo/reference_data/gencode.v44.basic.annotation_df.txt"
    "/directflow/SCCGGroupShare/projects/anncuo/reference_data/gencode.v44.basic.annotation_df.txt"
)
gene_info = pd.read_csv(gencode_text_file)

# Add gene info to AnnData
# set gene ids as indices
adata.var.index = [gene for gene in adata.var["gene_ids"]]
genes_cellranger = adata.var.index
# Remove gene version
gene_info.index = [gene.split(".")[0] for gene in gene_info["gene_id"]]
# Only consider genes in AnnData object
gene_info_df = gene_info[gene_info.index.isin(genes_cellranger)]
# Add info
adata.var = pd.concat([adata.var, gene_info_df], axis=1)

# remove batch from cell index (handling cases of a hyphen in the sample name)
adata.obs.index = [re.sub(r"-[0-9]+$", "", cell) for cell in adata.obs.index]

# perform QC filtering

# remove platelet, erythrocyte, and doublets as predicted by scPred
adata = adata[
    ~adata.obs["wg2_scpred_prediction"].isin(["Platelet", "Eryth", "Doublet"])
]

# Apply qc thresholds to remove low quality cells
print("number of cells pre-QC, pre-donor filtering:")
print(adata.obs.shape[1])
print(adata.obs.columns)

adata = adata[
    (adata.obs["pct_counts_mt"] < 20)
    & (adata.obs["n_genes_by_counts"] > 1000)
    & (adata.obs["n_genes_by_counts"] < 10000)
    & (adata.obs["total_counts"] > 800),
    :,
]


print("number of cells post-QC: before donor-filtering")
print(adata.obs.shape[1])

# get a subset of the metadata for UMAP QC plots
adata.obs.to_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_filtered_min1000genes_all_donors.csv"
)

# write
out_file = f"{out_dir}300_libraries/300_libraries_concatenated_filtered_all_donors.h5ad"
adata.write(out_file)

# DONOR-level QC:

# Remove any donors that have fewer than 100 cells after all previous QC is applied
# after applying this filter, 298 libraries remain. S0155 and S0156
donor_cell_counts = adata.obs["cpg_id"].value_counts()
donors_keep = donor_cell_counts.index[donor_cell_counts >= 100]
adata = adata[adata.obs["cpg_id"].isin(donors_keep)]

# remove data for donors with abnormal cell type composition
adata = adata[
    ~adata.obs["cpg_id"].isin(
        [  # (very high proportion of b intermediate)
            "CPG309724",
            "CPG310938",
            "CPG312025",
            "CPG315986",
            "CPG247973",
            "CPG249177",
            "CPG251793",
            # abnormal cell type distribution (other)
            "CPG252494",
            "CPG254169",
            "CPG254318",
            "CPG255760",
            "CPG249904",
        ]
    )
]

# remove data for donors that failed WGS QC
wgs_qc_fails = pd.read_csv(
    "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/saige-qtl_tenk10k-genome-2-3-eur_all_samples_to_drop.csv"
)

adata = adata[~adata.obs["cpg_id"].isin(wgs_qc_fails["s"])]

# Save the donor-filtered object (used in SAIGE eqtl analysis)

# after filtering on n_cells, 298 libraries remain.

# should be a copy of:
# "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_filtered_min1000genes.csv"
adata.obs.to_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/298_libraries_cell_metadata_filtered_min1000genes.csv"
)

print("number of cells post-QC: after donor-filtering")
print(adata.obs.shape[1])

# NOTE: This is same object as:
# 300_libraries/300_libraries_concatenated_filtered.h5ad
out_file = f"{out_dir}300_libraries/298_libraries_concatenated_filtered.h5ad"
adata.write(out_file)
