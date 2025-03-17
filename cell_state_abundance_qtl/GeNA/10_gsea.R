library(fgsea)
library(msigdbr)
library(tidyverse)
library(data.table)
library(glue)
library(patchwork)

celltype <- "Monocyte"
resolution <- "major_cell_types"
analysis_name <- "no_expr_pc_covars"

# Retrieve human H (hallmark) gene set
msigdbr_df_hallmark <- msigdbr(species = "human", category = "H")
msigdbr_df_immune <- msigdbr(species = "human", category = "C7")

msigdbr_df_immune %>%
    pull(gs_subcat) %>%
    unique()
# fixing format to work with fgsea
pathways_hallmark <- split(x = msigdbr_df_hallmark$ensembl_gene, f = msigdbr_df_hallmark$gs_name)
# pathways_immune <- split(x = msigdbr_df_immune$ensembl_gene, f = msigdbr_df_immune$gs_name)

# read in the gene correlations for csaQTLs

var_gene_correlations <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_vargene_corrs.csv")) %>%
    rename("gene_id" = V1)

csaQTL_npheno_corr_gene_ranks <- list()

for (col in var_gene_correlations %>%
    select(-gene_id) %>%
    colnames()) {
    rank_vec <- setNames(var_gene_correlations[[col]], var_gene_correlations[["gene_id"]]) %>% sort()
    csaQTL_npheno_corr_gene_ranks[[col]] <- rank_vec
}

# monocyte case examples
variant_list <- c("12:69350234:C:A", "12:9953308:T:TG")

for (variant in variant_list) {
    fgsea_hallmark <- fgsea(pathways = pathways_hallmark, csaQTL_npheno_corr_gene_ranks[[glue("npheno_{variant}")]]) %>%
        filter(padj <= 0.05) %>%
        arrange(padj)

    # save the results
    fgsea_hallmark %>%
        as_tibble() %>%
        mutate(leadingEdge = unlist(map(leadingEdge, ~ str_c(.x, collapse = ",")))) %>%
        write_tsv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/fgsea/{celltype}_npheno_{variant}_fgsea_sig.tsv"))

    plot_list <- list()

    for (pathway in fgsea_hallmark$pathway) {
        plot_list[[pathway]] <- plotEnrichment(
            pathways_hallmark[[pathway]],
            csaQTL_npheno_corr_gene_ranks[[glue("npheno_{variant}")]]
        ) +
            labs(title = pathway)
    }

    gsea_plots <- plot_list %>%
        wrap_plots(ncol = 3)

    gsea_plots %>% ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/gsea/{celltype}_{variant}_enrichplots.pdf"),
        width = 15, height = length(gsea_plots)
    )

    # save the genes contained in each gene set
    sig_gene_set_genes <- msigdbr_df_hallmark %>%
        filter(gs_name %in% fgsea_hallmark$pathway)
    sig_gene_set_genes %>%
        write_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/fgsea/{celltype}_npheno_{variant}_fgsea_genes.tsv"))
}
