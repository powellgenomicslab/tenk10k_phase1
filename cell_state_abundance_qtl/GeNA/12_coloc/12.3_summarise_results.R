library(tidyverse)
library(data.table)

# coloc_results_dir <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/"

# results_file_list <- list.files(, full.names=TRUE)

coloc_all_blood <- Sys.glob("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/*/blood_traits/*.tsv") %>%
    map(\(f)    {
        pheno_name <- f %>%
            str_remove("^.*/") %>%
            str_remove("_[0-9]+_[0-9]+\\.tsv")

        coloc_res <- fread(f) %>%
            mutate(gwas_pheno_name = pheno_name) %>%
            return()
    }) %>%
    rbindlist(use.names = TRUE)


coloc_all_blood %>%
    filter(PP.H4.abf > 0.8)

coloc_all_blood %>%
    pull(nsnps_coloc_tested) %>%
    max()


coloc_all_disease <- Sys.glob("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/*/disease_traits/*.tsv") %>%
    map(\(f)    {
        pheno_name <- f %>%
            str_remove("^.*/") %>%
            str_remove("_[0-9]+_[0-9]+\\.tsv")

        coloc_res <- fread(f) %>%
            mutate(gwas_pheno_name = pheno_name) %>%
            return()
    }) %>%
    rbindlist(use.names = TRUE)

coloc_all_disease %>%
    filter(PP.H4.abf > 0.8)

coloc_all_disease %>%
    arrange(desc(PP.H4.abf))


# test <- fread("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/coloc_results/gymrek-ukbb-eosinophil_count/CD4_Naive/gene_summary_result.csv")

# eqtl_coloc_results
