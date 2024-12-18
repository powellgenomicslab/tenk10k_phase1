# conda activate r_python_dev

library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(Seurat)

source("plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")


# 📗 Read in the data ----

# read in the cell-level metadata
# cell_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_cell_metadata_subset.csv") %>%
#     rename("barcode" = 1)

# cell_cycle_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/240_libraries_cellcyle_phase.csv") %>%
#     rename("barcode" = 1)

cell_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_cell_metadata_subset_filtered_reanalysed.csv") %>%
    rename("barcode" = 1)

cell_cycle_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/240_libraries_cellcyle_phase.csv") %>%
    rename("barcode" = 1)

# read in the UMAP coordinates
# umap_coords <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_harmony_umap_covs.csv") %>%
#     rename("barcode" = 1)

umap_coords <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_harmony_umap_coords_filtered_reanalysed.csv") %>%
    rename("barcode" = 1)

plot_data <- umap_coords %>%
    left_join(cell_metadata) %>%
    left_join(cell_cycle_metadata) %>%
    rename("UMAP 1" = UMAP1, "UMAP 2" = UMAP2)

# ⚙️ Functions ----

ggFeaturePlot <- function(
    data,
    feature,
    save = FALSE,
    path = NULL,
    size = 0.1,
    stroke = 0,
    ...) {
    data <- arrange(data, feature) # plot the lowest expressing cells first and highest expression last

    umap_ratio <- square.ratio(data[, c("UMAP 1", "UMAP 2")])
    feature_plot <- ggplot(data) +
        geom_point(
            mapping = aes(x = .data[["UMAP 1"]], y = .data[["UMAP 2"]], colour = .data[[feature]]),
            size = size, stroke = stroke, ...
        ) +
        labs(title = feature) +
        umap_ratio +
        theme_bw() +
        theme(
            legend.position = "right",
            legend.text.align = 0,
            legend.key.size = unit(5, "mm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            plot.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            ...
        ) +
        scale_color_gradientn(colours = c("lightgrey", "red"))

    if (save) {
        ggsave(feature_plot,
            filename = glue("{path}/ggfeature_{feature}.png"),
            width = 9, height = 5,
            dpi = 1300
        )
    } else {
        return(feature_plot)
    }
}


# 📊 make the plots ----

# UMAP colored by cell type
plot_data <- plot_data %>%
    mutate(wg2_scpred_prediction = factor(wg2_scpred_prediction, levels = tenk_color_pal$cell_type))

cohort_plot <- plot_data %>%
    ggUMAPplot(group.by = "cohort")

cohort_plot %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggumap_cohort.png",
        width = 9, height = 5,
        dpi = 1300,
    )


# Umap colored by ell cycle phase

plot_data %>%
    ggUMAPplot(group.by = "phase") %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggumap_cell_cycle_phase.png",
        width = 9, height = 5,
        dpi = 1300
    )

# Quality metric feature plots
# quality_metrics <- c("n_genes_by_counts", "total_counts", "pct_counts_mt", "S_score", "G2M_score")
quality_metrics <- c("n_genes_by_counts", "total_counts", "pct_counts_mt")

quality_metrics %>%
    purrr::walk(\(metric) ggFeaturePlot(
        data = plot_data,
        feature = metric,
        size = 0.1,
        save = TRUE,
        path = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/"
    ))


# Violin plots ----

# TODO: violin plot function

# ggViolin <- function(data) {

# }

# make violin plots for each cell type / qc metric

n_genes_by_counts_by_cell_types <- plot_data %>%
    pivot_longer(cols = {{ quality_metrics }}, names_to = "qc_metric", values_to = "qc_metric_value") %>%
    ggplot(aes(x = wg2_scpred_prediction, y = qc_metric_value, fill = wg2_scpred_prediction, color = wg2_scpred_prediction)) +
    facet_wrap(~qc_metric, scales = "free_y", ncol = 1) +
    geom_violin(width = 1) +
    scale_fill_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
    scale_color_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

n_genes_by_counts_by_cell_types %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggVln_n_genes_by_counts_by_cell_types.png",
        width = 14,
        height = 7,
        dpi = 300
    )

# TODO: make violin plots for each sequencing batch

# export ncells per sample

# cell_metadata %>%
#     group_by(cpg_id, cohort) %>%
#     summarise(total_barcodes = n()) %>%
#     write_csv("240_lib_all_individuals_ncells.csv")
