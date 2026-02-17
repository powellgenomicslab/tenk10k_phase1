import sys
import numpy as np
import pandas as pd
import cna
import multianndata as mad
import scanpy as sc

# celltype = sys.argv[1]
# resolution = sys.argv[2]
# analysis_name = sys.argv[3]

resolution = "major_cell_types"
analysis_name = "no_expr_pc_covars"
celltype = "Monocyte"

variant_list = "12:69350234:C:A,12:9953308:T:TG"  # add the variants here, the script will fetch gene sets that are correllated with the neighbohood level phenotype of these variants
variant_list = variant_list.split(",")

# read in the single cell data
madata = cna.read(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_scDataObject.dimreduc.pca.gena_pheno.h5ad"
)

for variant in variant_list:
    # read in and parse the top enriched gene sets for variant of interest
    fgsea = pd.read_csv(
        f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/fgsea/{celltype}_npheno_{variant}_fgsea_sig.tsv",
        sep="\t",
    )

    # Parse the fgsea results turn gene set back to list
    fgsea["leadingEdge"] = fgsea["leadingEdge"].str.split(",")

    # read in the file containing the genes for each geneset
    geneset_genes = pd.read_csv(
        f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/fgsea/{celltype}_npheno_{variant}_fgsea_genes.tsv",
    )

    # compute the module score for each geneset

    for gs in fgsea["pathway"].values:

        gl = geneset_genes.loc[geneset_genes["gs_name"] == gs, "ensembl_gene"].tolist()
        print(
            len(set(gl).intersection(set(madata.var.index)))
        )  # check how many genes in the gene set are also in the highly variable genes
        sc.tl.score_genes(madata, gene_list=gl, score_name=gs)

    madata.obs.columns
    # plot some UMAPs and see how well the neighborhood-level phenotypes align with the module scores

    sc.settings.figdir = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/gsea/"

    sc.pl.umap(
        madata,
        color=["wg2_scpred_prediction", f"npheno_{variant}"]
        + fgsea["pathway"].values.tolist(),
        vmin="p5",
        vmax="p95",
        save=f"_{celltype}_{variant}_umaps.png",
    )

madata.obs.to_csv(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_neighbourhood_pheno_modulescores.csv"
)
