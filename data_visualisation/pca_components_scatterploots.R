# look at correlation between PCA components and metadata + qc metrics (used in SAIGE-QTL model)
# useful to see if the PCA components are capturing the effects of QC metrics and potential confounders such as QC metrics, cohort, etc.

library(tidyverse)
library(glue)


celltype <- "CD4_TCM"
# pca_version = "" # pca with harmony + regress out %MT and total counts
# pca_version = "_no_harmony_with_regression"
pca_version = "_no_harmony_no_regression"
# data_type <- "filtered" # use filtered data, currently filtered with 600>n_genes>6000, <20% MT, >1000 total_counts

cell_metadata <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_subset_filtered_min1000genes.csv")) %>%
    rename("barcode" = 1)

celltype_pca_df <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/cpg_cell_covs_filtered/{celltype}_expression_pcs{pca_version}.csv")) %>%
    rename("barcode" = 1) %>%
    select(c("barcode", paste0("PC", 1:20)))
    
pca_df_long <- celltype_pca_df %>%
    pivot_longer(paste0("PC", 1:20), names_to = "PCA_component", values_to = "PCA_value")

# combine pca components df with cell metadata
plot_data <- pca_df_long %>%
    left_join(cell_metadata, by = "barcode")

# make scatter plots

make_scatter <- function(plot_data, metric, celltype, pca_version) {
    scatter_plots <- plot_data %>%
        ggplot(aes(x = PCA_value, y = !!sym(metric))) +
        facet_wrap(~PCA_component) +
        geom_point(
            size = 0.1,
            stroke = 0
        ) +
        theme_classic() +
        theme(aspect.ratio = 1)

    scatter_plots %>%
        ggsave(
            filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/pca_scatter_plots/{celltype}_300_libraries_{metric}_vs_PCs_scatterplots{pca_version}_20_pcs.png"),
            width = 15, height = 15
        )
}

metrics_to_plot = c("total_counts", "n_genes", "pct_counts_mt", "pct_counts_ribo", "cohort", "cellbender_background_fraction")

metrics_to_plot %>% map(\(metric) {
    make_scatter(plot_data = plot_data, metric = metric, celltype = celltype, pca_version = pca_version)
})

# make_scatter(plot_data = plot_data, metric = "total_counts", celltype = celltype)
# make_scatter(plot_data = plot_data, metric = "n_genes", celltype = celltype)
# make_scatter(plot_data = plot_data, metric = "n_genes", celltype = celltype)

# calculate correllation between each PC and each qc metric

# get_corr <- function() { }

metric <- "total_counts"

metric <- cell_metadata[cell_metadata$barcode == pca_df_long$barcode, metric]




cor(method="spearman")
