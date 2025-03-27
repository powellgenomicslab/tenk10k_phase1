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
csaQTL_npheno <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/major_cell_types/no_expr_pc_covars/{celltype}_neighbourhood_pheno.csv")) %>%
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
npheno_casexample <- ggFeaturePlot(plot_data, feature = "npheno_15:39687137:C:T") + labs(title = "15:39687137:C:T")

npheno_lyz <- ggFeaturePlot(plot_data, feature = "npheno_12:69350234:C:A") + labs(title = "12:69350234:C:A")
npheno_clec12a <- ggFeaturePlot(plot_data, feature = "npheno_12:9953308:T:TG") + labs(title = "12:9953308:T:TG")

# npheno_casexample %>%
#     ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_featureplot_casexample_npheno.png", width = 4, height = 4)

npheno_casexample %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_featureplot_casexample_npheno.pdf"), width = 4, height = 4)

# combined plots

combined <- umap_celltypes + npheno_casexample + plot_layout(ncol = 2)

combined %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_combined_umaps.pdf"), width = 8, height = 4)


# combined with lyz example
combined_lyz <- umap_celltypes + npheno_casexample + npheno_lyz + npheno_clec12a + plot_layout(ncol = 1)

combined_lyz %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_combined_umaps_lyz.pdf"), width = 4, height = 16)
