import scanpy as sc
import pandas as pd
import sys
import scanpy.external as sce

sc.settings.verbosity = 4  # show logging output


def run_umap(filename, n_pcs, n_neighors=15):
    """_summary_

    Args:
        filename (str):  file path to scanpy input h5ad file
        n_pcs (int): number of principal components to use in nearest neighbour graph construction
        n_neighbors(int): the number of neighbors to use (k) in the nearest neighbor graph
    """
    adata = sc.read(filename, cache=True)
    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighors=n_neighors)
    sc.tl.umap(adata)
    # save the UMAP coords for plotting in R
    pd.DataFrame(
        adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
    ).to_csv(
        f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/umap_coords_{n_pcs}_pcs.csv"
    )


def test_umap_hvg_selection(filename, n_hvg=3000, cohort=None):
    """test different hvg selection methods and see how it affects the UMAP

    Args:
        filename (str): path to input h5ad
        n_hvg (int): number of hvg to use
        cohort (str): if
    """

    adata = sc.read(filename, cache=True)
    adata = adata.raw.to_adata()

    valid_cohort = ["TOB", "BioHEART"]

    if cohort is None:
        out_umap_csv = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/umap_coords_30_pcs_{n_hvg}_hvg_outliers_removed.csv"

    elif cohort in valid_cohort:
        adata = adata[adata.obs["cohort"] == cohort]
        out_umap_csv = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/umap_coords_30_pcs_{n_hvg}_hvg_{cohort}_only_outtliers_removed.csv"

    else:
        raise ValueError(f"Invalid argument. Expected one of: {valid_cohort}")

    adata = adata[
        ~adata.obs["cpg_id"].isin(
            [
                "CPG309724",
                "CPG310938",
                "CPG312025",
                "CPG315986",
                "CPG247973",
                "CPG249177",
                "CPG251793",
            ]
        )
    ]

    # only consider genes with more than 3 counts across all cells
    sc.pp.filter_genes(adata, min_counts=3)

    # Use Anna's preprocessing steps from prepare_cell_covs.py
    # Note: I've changed the hvg selection method
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_hvg, batch_key="sequencing_library"
    )
    # adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(adata, max_value=10)

    # Dim Reduction
    sc.tl.pca(adata, svd_solver="arpack", n_comps=30)

    # perform Harmony integration to remove technical batch effects between pools
    sce.pp.harmony_integrate(adata, "sequencing_library")

    # UMAP on the Harmony principal components
    adata.obsm["X_pca"] = adata.obsm["X_pca_harmony"]
    sc.pp.neighbors(
        adata, n_pcs=30
    )  # NOTE: chose 30 PC's a based on the above analysis + elbow plot
    sc.tl.umap(adata)

    # export umap
    pd.DataFrame(
        adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index
    ).to_csv(out_umap_csv)


def main():
    # n_pcs_list = list(range(20, 46, 5))
    # n_pcs_list = list(range(10, 16, 5))
    # # Test UMAP over a range of PCs
    # for n_pcs in n_pcs_list:
    #     run_umap(
    #         filename="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony_leiden_filtered_reanalysed.h5ad",
    #         n_pcs=n_pcs,
    #     )

    # n_hvg = int(sys.argv[1])

    # Note: run test_umap_hvg_selection in parallel / different jobs because it takes a long time to run

    # test_umap_hvg_selection(
    #     filename="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony_leiden_filtered_reanalysed.h5ad",
    #     n_hvg=n_hvg,
    # )

    # test_umap_hvg_selection(
    #     filename="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony_leiden_filtered_reanalysed.h5ad",
    #     cohort="BioHEART",
    # )

    # test_umap_hvg_selection(
    #     filename="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony_leiden_filtered_reanalysed.h5ad",
    #     cohort="TOB",
    # )

    test_umap_hvg_selection(
        filename="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony_leiden_filtered_reanalysed.h5ad",
    )


if __name__ == "__main__":
    main()
