import sys
import pandas as pd
import multianndata as mad
import cna
import scanpy as sc
from matplotlib import pyplot as plt

celltype = sys.argv[1]
resolution = sys.argv[2]

outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

# ---------------------------
# sample-level pseudobulk PCA
# ---------------------------

# read in the data

adata = sc.read(
    f"{outdir}/data/h5/{resolution}/{celltype}_scanpy.h5ad",
    cache=True,
)

donor_mapping = adata.obs["cpg_id"]

adata.X = adata.layers["counts"]

# aggregate counts per-sample
pbcounts = adata.to_df().join(donor_mapping).groupby("cpg_id").agg("sum").sort_index()

# prepare the pseudobulk metadata
obs_pb = (
    adata.obs[["cpg_id", "cohort"]]
    .copy()
    .drop_duplicates()
    .set_index("cpg_id")
    .sort_index()
)

obs_pb["BioHEART"] = obs_pb["cohort"].map({"BioHEART": 1, "TOB": 0})
# obs_pb = obs_pb.drop("cohort", axis=1)

# make pseuodbulk anndata
adata_pb = sc.AnnData(X=pbcounts, obs=obs_pb, var=adata.var)

# re-calculate qc metrics
sc.pp.calculate_qc_metrics(
    adata_pb, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pp.normalize_total(adata_pb, target_sum=1e6)
sc.pp.log1p(adata_pb)
sc.pp.highly_variable_genes(adata_pb, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_pb = adata_pb[:, adata_pb.var.highly_variable]

# NOTE Leave out the step regressing out total_counts and pct_counts_mt
# I'm not including them as covars in the GeNA model, so hopefully they are captured

# UPDATE: actually I might want to include them in the model instead of the expression PCs
# im concerned expression PCs might also potentially capture cell type abundance shifts due to pseudobulking
# in which case, I would prefer not to include them as covariates, as these abundance shifts are what i'm trying
# to find

sc.pp.scale(adata_pb, max_value=10)
sc.tl.pca(adata_pb, svd_solver="arpack")

with plt.rc_context(
    {"figure.figsize": (5, 5)}
):  # Use this to set figure params like size and dpi
    sc.pl.pca(
        adata_pb,
        color=[
            "cohort",
            "cohort",
            "total_counts",
            "total_counts",
            "n_genes_by_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "pct_counts_mt",
        ],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3), (0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        show=False,
    )
    plt.savefig(
        f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/pca/{resolution}_{celltype}_pseudobulk_pcaplot.png",
        bbox_inches="tight",
    )

with plt.rc_context({"figure.figsize": (5, 5)}):
    sc.pl.pca_variance_ratio(adata_pb, n_pcs=50, log=True)
    plt.savefig(
        f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/pca/{resolution}_{celltype}_pseudobulk_elbow.png",
        bbox_inches="tight",
    )

# NOTE seems like typically sequencing depth and cohort are separated across PCA dim1 and dim2

donor_pca = pd.DataFrame(adata_pb.obsm["X_pca"])
donor_pca.index = adata_pb.obs.index
donor_pca.columns = [f"PC{i+1}" for i in donor_pca.columns]

donor_covariates = pd.DataFrame(
    adata_pb.obs[["BioHEART", "total_counts", "pct_counts_mt"]]
)


# donor_covariates["total_counts"].corr(donor_covariates["pct_counts_mt"])


# ----
# add donor-PCA to the madata sample matrix
# ----

madata = cna.read(
    f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.h5ad"
)

print(madata.samplem)

madata.samplem = madata.samplem.join(donor_pca).join(donor_covariates)

print(madata.samplem)

print(madata)

# over-write output from previous step
madata.write(f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.pca.h5ad")
