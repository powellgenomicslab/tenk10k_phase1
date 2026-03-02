import pandas as pd
import click
import os
import re
import sys
import logging
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

# 1) get lead csaQTL GWAS loci
# 2) Extract the eQTL data for the cis-window (lead csaQTL GWAS loci +/- window; default: 100kB)
# 3) Run coloc for each lead csaQTL (For significant egenes pval > 0.05)
# 4) Write the results to a TSV file


def coloc_runner(
    eqtl_results_cis_gene,
    csaqtl_file_path,
    celltype_csaqtl,
    celltype_eqtl,
    outdir,
    common_maf_threshold,
    gene,
    save=True,
):
    """Modified function from: https://github.com/populationgenomics/saige-tenk10k/blob/main/coloc/coloc_runner.py"""

    ro.r("suppressMessages(library(coloc))")
    ro.r("suppressMessages(library(tidyverse))")

    eqtl = eqtl_results_cis_gene
    # eqtl['chr'] = eqtl['CHR'].apply(lambda chr: f'chr{chr}')
    eqtl["beta"] = eqtl["BETA"]
    eqtl["se"] = eqtl["SE"]
    # eqtl['position'] = eqtl['POS']
    eqtl["snp"] = eqtl["MarkerID"].apply(lambda snp: "chr" + snp.replace(":", "_"))
    # while I figure out if it's easy to extract sdY, give N and MAF instead
    # https://chr1swallace.github.io/coloc/articles/a02_data.html#what-if-i-dont-have-sdy
    eqtl["MAF"] = eqtl["AF_Allele2"].apply(lambda af: min(af, (1 - af)))
    # subset to only results for variants above chosen MAF
    eqtl = eqtl[eqtl["MAF"] > common_maf_threshold]
    with (ro.default_converter + pandas2ri.converter).context():
        eqtl_r = ro.conversion.get_conversion().py2rpy(eqtl)
    ro.globalenv["eqtl_r"] = eqtl_r
    ro.globalenv["gene"] = gene
    ro.r(
        """
    eqtl_r = eqtl_r %>% filter(!is.na(beta))
    eqtl_r = eqtl_r %>% distinct(snp, .keep_all = TRUE)
    eqtl_r$varbeta = eqtl_r$se**2
    eqtl_r$position = eqtl_r$pos
    eqtl_r = eqtl_r %>% select(beta, varbeta, position, snp, N, MAF)

    eqtl_r = eqtl_r %>% as.list()
    eqtl_r$type = 'quant'

    """,
    )

    csaqtl = pd.read_csv(
        csaqtl_file_path,
        sep="\t",
    )
    csaqtl["chr"] = csaqtl["#CHROM"].apply(lambda chr: f"chr{chr}")
    csaqtl["beta"] = csaqtl["BETA"]
    csaqtl["se"] = csaqtl["SE"]
    csaqtl["position"] = csaqtl["POS"]
    csaqtl["snp"] = csaqtl["ID"].apply(lambda snp: "chr" + snp.replace(":", "_"))
    # while I figure out if it's easy to extract sdY, give N and MAF instead
    # https://chr1swallace.github.io/coloc/articles/a02_data.html#what-if-i-dont-have-sdy
    csaqtl["MAF"] = csaqtl["A1_FREQ"]
    csaqtl["N"] = csaqtl["OBS_CT"]  # n_samples
    # subset to only results for variants above chosen MAF (0.05 by default)
    csaqtl = csaqtl[csaqtl["MAF"] > common_maf_threshold]
    with (ro.default_converter + pandas2ri.converter).context():
        csaqtl_r = ro.conversion.get_conversion().py2rpy(csaqtl)

    lead_snp = "_".join(
        csaqtl_file_path.split("/")[-1].split("_")[2:4]
    )  # Parse lead SNP from the

    ro.globalenv["csaqtl_r"] = csaqtl_r
    ro.globalenv["lead_snp"] = lead_snp
    ro.r(
        """
    csaqtl_r = csaqtl_r %>% filter(!is.na(beta))
    csaqtl_r = csaqtl_r %>% distinct(snp, .keep_all = TRUE)
    csaqtl_r$varbeta = csaqtl_r$se**2
    csaqtl_r$position = csaqtl_r$pos
    csaqtl_r = csaqtl_r %>% select(beta, varbeta, position, snp, N, MAF)

    csaqtl_r = csaqtl_r %>% as.list()
    csaqtl_r$type = 'quant'

    my.res <- coloc.abf(dataset1=eqtl_r,
                    dataset2=csaqtl_r)
    p_df <- data.frame(lead_snp, my.res$summary[1], my.res$summary[2], my.res$summary[3], my.res$summary[4], my.res$summary[5], my.res$summary[6])
    names(p_df) <- c('lead_snp', 'nsnps_coloc_tested','PP.H0.abf','PP.H1.abf','PP.H2.abf','PP.H3.abf','PP.H4.abf')
    """,
    )

    # convert to pandas df
    with (ro.default_converter + pandas2ri.converter).context():
        pd_p4_df = ro.conversion.get_conversion().rpy2py(ro.r("p_df"))

    # add cell type and chrom annotation to df
    pd_p4_df["celltype_csaqtl"] = celltype_csaqtl
    pd_p4_df["celltype_eqtl"] = celltype_eqtl
    pd_p4_df["chrom"] = csaqtl["chr"].iloc[0]

    if save:
        coloc_results_file = (
            f"{outdir}/{celltype_csaqtl}_{celltype_eqtl}_{lead_snp}_{gene}.tsv"
        )
        # write to GCS
        pd_p4_df.to_csv(
            coloc_results_file,
            sep="\t",
            index=False,
        )
    else:
        return pd_p4_df


@click.option(
    "--egenes-files-path",
    help="Path to the directory containing the eQTL eGenes files",
)
@click.option(
    "--celltype-eqtl",
    help="name of the minor cell type of csaQTLs being colocalised",
)
@click.option(
    "--celltype-csaqtl",
    help="name of the major cell type of csaQTLs being colocalised",
)
@click.option(
    "--csaqtl-file-path",
    help="Path to the csaQTL fixed-phenotype summary statistics files",
)
@click.option(
    "--cis-window",
    help="Window around lead csaQTL variant to test for colocalisations",
    default=100000,
)
@click.option(
    "--fdr-threshold",
    help="FDR threshold for eQTL ACAT p-values to select genes for colocalisation",
    default=0.05,
)
@click.option(
    "--common-maf-threshold",
    help="eQTL SNPs common MAF threshold",
    default=0.05,
)
@click.option("--outdir", help="Output directory")
@click.command()
def main(
    egenes_files_path: str,
    celltype_csaqtl: str,
    celltype_eqtl: str,
    csaqtl_file_path: str,
    common_maf_threshold: float,
    fdr_threshold: float,
    cis_window: int,
    outdir: str,
):
    os.makedirs(outdir, exist_ok=True)
    chr, pos = [
        int(i) for i in csaqtl_file_path.split("/")[-1].split("_")[2:4]
    ]  # parse lead SNP chr-pos from the file name
    lead_snp = f"{chr}_{pos}"
    # read in gene annotation file
    # var_table = pd.read_csv(gene_info_file)

    # read in eGenes file)
    egenes_file = f"{egenes_files_path}/{celltype_eqtl}/{celltype_eqtl}_all_cis_cv_gene_level_results.tsv"
    egenes_results = pd.read_csv(
        egenes_file,
        sep="\t",
    )
    egenes_results["chr"] = egenes_results["top_MarkerID"].apply(
        lambda snp: "chr" + snp.split(":")[0]
    )
    egenes_results = egenes_results[
        egenes_results["ACAT_p"] < fdr_threshold
    ]  # filter for sc-eQTLs with p-value < fdr_threshold
    # NOTE shouldn't this be qvalue not pvalue?

    # now read in the full eqtl results for all SNPs
    # TODO check with Anna - do these have all SNPs or only cis SNPs in 100kb cis window
    eqtl_results_all_snps = pd.read_csv(
        f"{egenes_files_path}/{celltype_eqtl}/{celltype_eqtl}_common_all_cis_raw_pvalues.tsv",
        sep="\t",
        index_col=0,
    ).rename(columns={"CHR": "chr", "POS": "position"})
    eqtl_results_all_snps["chr"] = [f"chr{chr}" for chr in eqtl_results_all_snps["chr"]]

    # Filter eQTL results to include only SNPs within cis region around csaQTL lead-snp and only significant egenes based on the above egene table
    eqtl_results_cis_egenes = eqtl_results_all_snps.loc[
        (eqtl_results_all_snps["chr"] == f"chr{chr}")
        & (eqtl_results_all_snps["position"] > pos - cis_window)
        & (eqtl_results_all_snps["position"] < pos + cis_window)
        & (eqtl_results_all_snps["gene"].isin(egenes_results["gene"])),
        :,
    ]

    if len(eqtl_results_cis_egenes["gene"].unique()) == 0:
        logging.info(
            "No eQTL data for " + lead_snp + " in the cis-window: skipping...."
        )

    # loop through genes with < 1 significant eQTL SNP in the cis-window around the lead csaQTL SNP
    for gene in eqtl_results_cis_egenes["gene"].unique():

        # filter eQTL summstats to only include results for one gene
        eqtl_results_cis_gene = eqtl_results_cis_egenes.loc[
            eqtl_results_cis_egenes["gene"] == gene, :
        ]  # pass this df to coloc function instead of saving
        # celltype_eqtl=f'{outdir}/eqtl_cis_snps/{celltype_eqtl}_{lead_snp}_{gene}.tsv'

        coloc_results_file = (
            f"{outdir}/{celltype_csaqtl}_{celltype_eqtl}_{lead_snp}_{gene}.tsv"
        )

        if os.path.exists(coloc_results_file):
            logging.info(
                f"Output file for {gene} and {lead_snp} for {celltype_csaqtl}/{celltype_eqtl} already exists: skipping...."
            )

        else:
            # run coloc
            coloc_runner(
                eqtl_results_cis_gene=eqtl_results_cis_gene,
                celltype_eqtl=celltype_eqtl,
                csaqtl_file_path=csaqtl_file_path,
                celltype_csaqtl=celltype_csaqtl,
                common_maf_threshold=common_maf_threshold,
                outdir=outdir,
                gene=gene,
            )


if __name__ == "__main__":
    main()
