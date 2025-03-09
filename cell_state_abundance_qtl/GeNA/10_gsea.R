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
pathways_immune <- split(x = msigdbr_df_immune$ensembl_gene, f = msigdbr_df_immune$gs_name)

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

# run fgsea enrichment
fgsea_hallmark <- fgsea(pathways = pathways_hallmark, csaQTL_npheno_corr_gene_ranks[["npheno_19:16331208:G:A"]]) %>%
    filter(padj <= 0.05) %>%
    arrange(padj)
# fgsea_immune <- fgsea(pathways = pathways_immune, csaQTL_npheno_corr_gene_ranks[["npheno_19:16331208:G:A"]]) %>%
#     filter(padj <= 0.05) %>%
#     arrange(padj)

fgsea_hallmark %>%
    mutate(leadingEdge = as.character(leadingEdge)) %>%
    write_tsv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/fgsea/{celltype}_npheno_19:16331208:G:A_fgsea_sig.tsv"))

HALLMARK_TNFA_SIGNALING_VIA_NFKB_enrichplot <- plotEnrichment(
    pathways_hallmark[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
    csaQTL_npheno_corr_gene_ranks[["npheno_19:16331208:G:A"]]
) +
    labs(title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

HALLMARK_INFLAMMATORY_RESPONSE_enrichplot <- plotEnrichment(
    pathways_hallmark[["HALLMARK_INFLAMMATORY_RESPONSE"]],
    csaQTL_npheno_corr_gene_ranks[["npheno_19:16331208:G:A"]]
) +
    labs(title = "HALLMARK_INFLAMMATORY_RESPONSE")

HALLMARK_IL2_STAT5_SIGNALING_enrichplot <- plotEnrichment(
    pathways_hallmark[["HALLMARK_IL2_STAT5_SIGNALING"]],
    csaQTL_npheno_corr_gene_ranks[["npheno_19:16331208:G:A"]]
) +
    labs(title = "HALLMARK_IL2_STAT5_SIGNALING")

gsea_plots <- HALLMARK_TNFA_SIGNALING_VIA_NFKB_enrichplot +
    HALLMARK_INFLAMMATORY_RESPONSE_enrichplot +
    HALLMARK_IL2_STAT5_SIGNALING_enrichplot + plot_layout(ncol = 1)

gsea_plots %>% ggsave(
    filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/gsea/{celltype}_KLF2_enrichplots.pdf"),
    width = 5, height = 6
)


# fgsea_immune <- fgsea(pathways = pathways_immune, csaQTL_npheno_corr_gene_ranks[["npheno_19:16331208:G:A"]]) %>%
#     filter(padj <= 0.05) %>%
#     arrange(padj)

# fgsea_immune
