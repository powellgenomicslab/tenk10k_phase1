import sys
import pandas as pd
import multianndata as mad
import cna

celltype = sys.argv[1]


outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

# ----
# Calculate the per-donor PCs (pseudobulk PC's)
# ----

celltype_pca_sc = pd.read_csv(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/cpg_cell_covs_filtered/{celltype}_expression_pcs_no_harmony_with_regression.csv",
    index_col=0,
)
metadata = pd.read_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/298_libraries_cell_metadata_filtered_min1000genes.csv",
    index_col=0,
)
barcode_donor_map = metadata.loc[
    metadata["wg2_scpred_prediction"] == celltype, ["cpg_id", "cohort"]
]

donor_mean_pca = (
    celltype_pca_sc.join(barcode_donor_map, how="left")
    .groupby(["cpg_id", "cohort"])
    .mean()
)

# rename the PC colnames
donor_mean_pca.columns = [f"mean_expr_{pc}" for pc in donor_mean_pca.columns]

donor_mean_pca = donor_mean_pca.reset_index(level="cohort")
donor_mean_pca["BioHEART"] = donor_mean_pca["cohort"].map({"BioHEART": 1, "TOB": 0})
donor_mean_pca = donor_mean_pca.drop(["cohort"], axis=1)

# ----
# add donor-PCA to the madata sample matrix
# ----

madata = cna.read(f"{outdir}/data/h5/{celltype}_scDataObject.dimreduc.h5ad")

print(madata.samplem)

madata.samplem = madata.samplem.join(donor_mean_pca)

print(madata.samplem)

print(madata)

# over-write output from previous step
madata.write(f"{outdir}/data/h5/{celltype}_scDataObject.dimreduc.h5ad")
