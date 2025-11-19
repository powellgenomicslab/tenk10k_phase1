library(tidyverse)
library(viridis)
library(scales)
library(glue)
library(scattermore)
library(ggbeeswarm)
library(patchwork)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/plot_functions.R")
celltype <- "Monocyte"
resolution <- "major_cell_types"

# 📚 read in the data ----

# UMAP data
umap <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/umaps/{resolution}/{celltype}_umap.csv")) %>%
    rename("barcode" = `...1`, "UMAP 1" = UMAP1, "UMAP 2" = UMAP2)
# neighborhood-level phenotype
csaQTL_npheno <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/major_cell_types/no_expr_pc_covars/{celltype}_neighbourhood_pheno_modulescores.csv")) %>%
    rename("barcode" = `...1`)
sc_metadata <- get_latest_metadata() %>%
    select(barcode, wg2_scpred_prediction)

sc_metadata <- sc_metadata %>%
    mutate(cell_type = factor(wg2_scpred_prediction %>% str_replace("_", " "), levels = tenk_color_pal$cell_type))

plot_data <- umap %>%
    left_join(csaQTL_npheno, by = "barcode") %>%
    left_join(sc_metadata, by = "barcode")

# plot UMAP colored by cell types ----

umap_celltypes <- ggUMAPplot(plot_data, group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type))
umap_celltypes %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_umap_celltypes.pdf"), width = 4, height = 4)

# plot neighbourhood-level phenotype over UMAP ---

# this is the example altering discrete populations
# npheno_casexample <- ggFeaturePlot(plot_data, feature = "npheno_15:39687137:C:T") + labs(title = "15:39687137:C:T")

# npheno_casexample %>%
#     ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_featureplot_casexample_npheno.png", width = 4, height = 4)

# npheno_casexample %>%
# ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_featureplot_casexample_npheno.pdf"), width = 4, height = 4)

# combined plots

# combined <- umap_celltypes + npheno_casexample + plot_layout(ncol = 2)

# combined %>%
#     ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_combined_umaps.pdf"), width = 8, height = 4)

# combined with lyz + CLEC12A examples

# these are some examples that have altered abundance of less discrete cell state phenotypes ----

# plot neighborhood-level phenotypes
npheno_lyz <- ggFeaturePlot(plot_data, feature = "npheno_12:69350234:C:A") + labs(title = "12:69350234:C:A")
npheno_clec12a <- ggFeaturePlot(plot_data, feature = "npheno_12:9953308:T:TG") + labs(title = "12:9953308:T:TG")

# combined_lyz <- umap_celltypes + npheno_casexample + npheno_lyz + npheno_clec12a + plot_layout(ncol = 1)

# combined_lyz %>%
#     ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_combined_umaps_lyz.pdf"), width = 4, height = 16)

# combined case examples for continuous cell state  phenotypes ----

# plot some of the top significant enriched gene sets for these phenotypes
module_score_TNFA <- ggFeaturePlot(plot_data, feature = "HALLMARK_TNFA_SIGNALING_VIA_NFKB") #+ labs(title = "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
module_score_INFLAMMATORY <- ggFeaturePlot(plot_data, feature = "HALLMARK_INFLAMMATORY_RESPONSE") #+ labs(title = "HALLMARK_INFLAMMATORY_RESPONSE")
module_score_IL2 <- ggFeaturePlot(plot_data, feature = "HALLMARK_IL2_STAT5_SIGNALING") #+ labs(title = "HALLMARK_IL2_STAT5_SIGNALING")
module_score_IL6 <- ggFeaturePlot(plot_data, feature = "HALLMARK_IL6_JAK_STAT3_SIGNALING") #+ labs(title = "HALLMARK_IL6_JAK_STAT3_SIGNALING")


# CLEC12A combined plot
combined_with_gsea_clec12a <- umap_celltypes + npheno_clec12a + module_score_TNFA + module_score_INFLAMMATORY + module_score_IL2 + module_score_IL6 + plot_layout(ncol = 3)
combined_with_gsea_clec12a %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_combined_umaps_clec12a_modulescores.pdf"), width = 12, height = 8)

# LYZ combined plot
combined_with_gsea_lyz <- umap_celltypes + npheno_lyz + module_score_TNFA + module_score_INFLAMMATORY + plot_layout(ncol = 3)
combined_with_gsea_lyz %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_combined_umaps_lyz_modulescores.pdf"), width = 12, height = 8)

# All module score plots together
combined_module_score_umaps <- module_score_TNFA + module_score_INFLAMMATORY + module_score_IL2 + module_score_IL6 + plot_layout(ncol = 2)
combined_module_score_umaps %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_combined_umaps_all_modulescores.pdf"), width = 6, height = 6)
