# mm activate r-python-dev

library(tidyverse)
library(data.table)
library(glue)

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"
celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt")

for (celltype in celltypes[!celltypes %in% c("Other", "ALL")]) {
    # 1. Get summary stats for the SNPs within cis window around each lead csaQTL SNP ----

    #' read in csaQTL summary statistics
    read_summstats <- function(celltype, file_path, analysis_name, resolution) {
        sumstats <- fread(file_path)
        sumstats[, `:=`(
            P = fifelse(as.numeric(P) == 0, .Machine$double.xmin * .Machine$double.eps, as.numeric(P))
            # celltype = celltype
        )]
        return(sumstats)
    }

    # read in csaQTL summary stats for all snps and for lead snps
    sumstats_all_snps <- read_summstats(celltype = celltype, file_path = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats.txt"), analysis_name = analysis_name, resolution = resolution)
    sumstats_lead_snps <- read_summstats(celltype = celltype, file_path = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats_lead_snps_MAF_0.05.tsv"), analysis_name = analysis_name, resolution = resolution)

    lead_snps <- sumstats_lead_snps$ID

    # save list of lead SNPs to iterate through later
    lead_snps %>%
        write_lines(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats_lead_snps_MAF_0.05.snps"))

    #' get the summary statistics for snps within cis window around the lead snps, save each to a separate file
    get_cis_summstats <- function(lead_snp, cis_window = 100000, out_dir) {
        chr <- sumstats_lead_snps[ID == lead_snp, `#CHROM`]
        pos <- sumstats_lead_snps[ID == lead_snp, POS]
        cis_snps <- sumstats_all_snps[`#CHROM` == chr & POS > pos - cis_window & POS < pos + cis_window]
        # cis_snps[,celltype := NULL]
        dir.create(out_dir, showWarnings = FALSE)
        output_file_tsv <- glue("{out_dir}/GeNA_sumstats_{chr}_{pos}_cis_snps.tsv")
        output_file_snps <- glue("{out_dir}/GeNA_sumstats_{chr}_{pos}_cis_snps.snps")
        fwrite(cis_snps, output_file_tsv, sep = "\t")
        fwrite(cis_snps[, .(ID)], output_file_snps, sep = "\t", col.names = FALSE) # create list of snps to filter
        return(glue("{chr}_{pos}"))
    }

    lead_snps_chr_pos <- lead_snps %>%
        map(\(lead_snp)        {
            get_cis_summstats(
                lead_snp = lead_snp, cis_window = 100000, out_dir = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/lead_snps_ciswindow/")
            )
        })

    # using chr_pos to name the files // iterate through later, because some of the indels are too long to save as file names.
    if (sum(duplicated(lead_snps_chr_pos)) > 0) {
        warning("Lead SNPs with same CHR-POS! Rename these files using snp ID's!")
    }
    lead_snps_chr_pos %>%
        write_lines(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats_lead_snps_MAF_0.05.snps"))

    # 2. reformat sample-level phenotype files for input into plink ----
    spheno <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_sample_pheno.csv"))
    # reformat ids to work with plink, remove other metadata
    spheno <- spheno[, .SD, .SDcols = patterns("id|spheno")][
        ,
        `:=`(
            `#FID` = 0,
            IID = id
        )
    ][
        ,
        id := NULL
    ]
    setcolorder(spheno, c(
        "#FID",
        "IID"
    ))

    # rename the sample phenotypes with chromosome/pos (remove ref:alt), as otherwise some are too long to be used in file names downstream
    new_col_names <- spheno %>%
        colnames() %>%
        str_match("^.*:[0-9]+|#FID|IID") %>%
        str_replace(":", "_")

    setnames(spheno, new_col_names)
    fwrite(spheno, glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_sample_pheno_for_plink.tsv"), sep = "\t", na = "NA", quote = FALSE)
}


# TODO:
# create a genotype file for the cis-snps
# see:  https://github.com/immunogenomics/GeNA-applied/blob/main/run_gwas/05.munge_geno_cissnps.sh


# check is A1 allele always the alternative allele? (NO it's the minor allele!)

# compare allele frequencies with alt_freqs file
