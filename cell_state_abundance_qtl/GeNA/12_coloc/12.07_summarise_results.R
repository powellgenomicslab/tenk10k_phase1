library(tidyverse)
library(data.table)
library(ggplot2)
library(rtracklayer)
# coloc_results_dir <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/"

# results_file_list <- list.files(, full.names=TRUE)
# ----
# summarise coloc results with GWAS
# ----

coloc_all_blood <- Sys.glob("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/*/blood_traits/*.tsv") %>%
    map(\(f)    {
        pheno_name <- f %>%
            str_remove("^.*/") %>%
            str_remove("_[0-9]+_[0-9]+\\.tsv")

        coloc_res <- fread(f) %>%
            mutate(gwas_pheno_name = pheno_name) %>%
            mutate(gwas_category = "blood_trait_gymrek") %>%
            return()
    }) %>%
    rbindlist(use.names = TRUE)

coloc_all_disease <- Sys.glob("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/*/disease_traits/*.tsv") %>%
    map(\(f)    {
        pheno_name <- f %>%
            str_remove("^.*/") %>%
            str_remove("_[0-9]+_[0-9]+\\.tsv")

        coloc_res <- fread(f) %>%
            mutate(gwas_pheno_name = pheno_name) %>%
            mutate(gwas_category = "disease_trait") %>%
            return()
    }) %>%
    rbindlist(use.names = TRUE)


coloc_all_vukovic <- Sys.glob("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/*/extra_traits/*.tsv") %>%
    map(\(f)    {
        pheno_name <- f %>%
            str_remove("^.*/") %>%
            str_remove("_[0-9]+_[0-9]+\\.tsv")

        coloc_res <- fread(f) %>%
            mutate(gwas_pheno_name = pheno_name) %>%
            mutate(gwas_category = "blood_trait_vukovic") %>%
            return()
    }) %>%
    rbindlist(use.names = TRUE)

coloc_all <- coloc_all_blood %>%
    bind_rows(coloc_all_disease) %>%
    bind_rows(coloc_all_vukovic)

coloc_sig_table <- coloc_all %>%
    filter(PP.H4.abf > 0.8)
# %>%
# pull(lead_snp) %>%
# unique()

# coloc_all$lead_snp %>% unique() %>% length()
# Exploratory Plots

plot_data <- coloc_all %>%
    pivot_longer(PP.H0.abf:PP.H4.abf, names_to = "Hypothesis", values_to = "Posterior probability")

barplot_summary <- plot_data %>%
    ggplot(aes(x = `Posterior probability`)) +
    # facet_grid(rows = vars(celltype), cols = vars(Hypothesis)) +
    facet_wrap(~Hypothesis, ncol = 5) +
    geom_bar(stat = "bin", bins = 10) +
    theme_bw()

barplot_summary %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/coloc/coloc_summary_bar.png", width = 10, height = 2.5)

barplot_by_ct <- plot_data %>%
    ggplot(aes(x = `Posterior probability`)) +
    facet_grid(rows = vars(celltype), cols = vars(Hypothesis)) +
    # facet_wrap(~Hypothesis, ncol = 5) +
    geom_bar(stat = "bin", bins = 10) +
    theme_bw()

barplot_by_ct %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/coloc/coloc_bar_by_ct.png", width = 10, height = 14)

barplot_by_pheno <- plot_data %>%
    ggplot(aes(x = `Posterior probability`)) +
    facet_grid(rows = vars(gwas_pheno_name), cols = vars(Hypothesis)) +
    # facet_wrap(~Hypothesis, ncol = 5) +
    geom_bar(stat = "bin", bins = 10) +
    theme_bw()

barplot_by_pheno %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/coloc/coloc_bar_by_pheno.png", width = 10, height = 20)

per_pheno_summary <- plot_data %>%
    filter(`Posterior probability` > 0.8) %>%
    ggplot(aes(x = gwas_pheno_name, fill = Hypothesis)) +
    facet_wrap(~Hypothesis, ncol = 1, scales = "free_y") +
    facet_grid(rows = vars(Hypothesis), cols = vars(gwas_category), scales = "free") +
    geom_bar(stat = "count") +
    labs(y = "N loci (PP > 0.8)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

per_pheno_summary %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/coloc/coloc_bar_by_pheno.png", width = 14, height = 10)

# ----
# summarise coloc results with eQTLs
# ----

coloc_all_eqtl <- Sys.glob("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/*/eqtl/*.tsv") %>%
    map(\(f) {
        ensg <- f %>%
            str_remove(".*_") %>%
            str_remove(".tsv")

        coloc_res <- fread(f) %>%
            mutate(ensg = ensg) %>%
            return()
    }) %>%
    rbindlist(use.names = TRUE)

coloc_sig_eqtl <- coloc_all_eqtl %>%
    filter(PP.H4.abf > 0.8)

# get the gene names
gencode_path <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/reference_data/gencode.v44.basic.annotation.gtf.gz"
gencode <- rtracklayer::import(gencode_path)
gencode_genes <- gencode[gencode$type == "gene", c("gene_id", "gene_name")] # filter down to just genes
gencode_genes$gene_id_noversion <- gencode_genes$gene_id %>% str_remove("\\.[0-9]+") # remove version numbers to match eqtl results
# intersect(gencode_genes$gene_id_noversion, all_eqtl_fdr_sig_gr$gene)
gene_name_map <- gencode_genes %>%
    as.data.frame() %>%
    select(gene_id_noversion, gene_name)

coloc_sig_eqtl <- coloc_sig_eqtl %>%
    left_join(gene_name_map, by = c("ensg" = "gene_id_noversion"))

coloc_sig_eqtl %>%
    fwrite("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/colo_sig_eqtl_results.tsv", sep = "\t")

# coloc_sig_eqtl$lead_snp %>% unique() %>% length()
# coloc_sig_eqtl$celltype_csaqtl %>% table()


# coloc_sig_eqtl

# GWAS is trait 1
# csaQTL is trait 2

############################################################################################################################################################################################################

# test <- fread("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/coloc_results/gymrek-ukbb-eosinophil_count/CD4_Naive/gene_summary_result.csv")
# test$nsnps_coloc_tested %>% max()
# # eqtl_coloc_results
# coloc_all_blood$nsnps_coloc_tested %>% max()
