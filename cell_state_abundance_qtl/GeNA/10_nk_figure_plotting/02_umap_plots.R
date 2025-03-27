library(tidyverse)
library(viridis)
library(scales)
library(glue)
library(scattermore)
library(ggbeeswarm)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")


# 📚 read in the data ----

# NK cells
# UMAP data
nk_umap <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/NK_cellstate/NK_cellstate_umap_subtypes.csv") %>%
    rename("barcode" = `...1`, "UMAP 1" = UMAP1, "UMAP 2" = UMAP2) %>%
    mutate(
        celltype = case_match(
            celltype,
            "NK1" ~ "NK_mature",
            "NK3" ~ "NK_adaptive",
            "NKint" ~ "NK_intermediate",
            .default = celltype
        )
    )
nk_csaQTL_npheno <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/major_cell_types/no_expr_pc_covars/NK_neighbourhood_pheno.csv") %>%
    rename("barcode" = `...1`)

plot_data <- nk_umap %>%
    left_join(nk_csaQTL_npheno, by = "barcode")

# plot UMAP colored by NK cell types ----

umap_celltypes <- ggUMAPplot(plot_data, group.by = "celltype")
umap_celltypes %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/NK_umap_celltypes.pdf", width = 4, height = 4)

# plot neighbourhood-level phenotype over UMAP ---

npheno_klf2 <- ggFeaturePlot(plot_data, feature = "npheno_19:16331208:G:A")

# npheno_klf2 %>%
#     ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/NK_featureplot_klf2_npheno.png", width = 4, height = 4)

npheno_klf2 %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/NK_featureplot_klf2_npheno.pdf", width = 4, height = 4)
