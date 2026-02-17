library(tidyverse)
library(data.table)
library(glue)

# reformats summary stats (.ma format) to work with the coloc pipeline
reformat_gwas <- function(pheno_name, filepath, outdir) {
    gwas_sumstats <- fread(filepath)

    gwas_formatted_for_coloc <- gwas_sumstats %>%
        mutate(
            chromosome = SNP %>% str_extract("[0-9]+:") %>% str_remove(":"),
            position = SNP %>%
                str_extract(":[0-9]+:") %>%
                str_remove_all(":"),
            snp = SNP %>% str_replace_all(":", "_")
        ) %>%
        mutate(
            chromosome = paste0("chr", chromosome),
            snp = paste0("chr", snp),
            varbeta = se**2,
        ) %>%
        select(
            chromosome,
            position,
            varbeta,
            beta = b,
            snp,
            MAF = freq,
            p_value = p
        )

    fwrite(gwas_formatted_for_coloc, glue("{outdir}/{pheno_name}.ma"), sep = "\t", na = "NA", quote = FALSE)
}

# run the function to save re-formatted summary statistics which can be used in coloc analysis
files <- Sys.glob("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/from_smr/*.ma")

pheno_names <- files %>%
    basename() %>%
    str_remove(".ma")


for (i in seq_along(files)) {
    reformat_gwas(pheno_names[i], files[i], "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/from_smr/reformatted_for_coloc")
}

# NOTES;

# freq is MAF. summstats already filtered to MAF > 0.01
# gwas_formatted_for_coloc$freq %>%
#     min()

# here is an example of how the formatting is for the other GWAS used downstream (preprocessed by Hope)
# gwas_formatted <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/white_british_alanine_aminotransferase_snp_str_gwas_results_hg38_chr1.tab.gz")
