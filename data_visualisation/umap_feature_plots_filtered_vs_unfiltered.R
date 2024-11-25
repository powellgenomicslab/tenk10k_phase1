# conda activate r_python_dev

library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(patchwork)

source("plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

data_type <- "unfiltered"

# 📗 Read in the data ----

# read in the cell-level metadata
# cell_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_cell_metadata_subset.csv") %>%
#     rename("barcode" = 1)

# cell_cycle_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/240_libraries_cellcyle_phase.csv") %>%
#     rename("barcode" = 1)

cell_metadata <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_subset_{data_type}.csv")) %>%
    rename("barcode" = 1)

cell_metadata <- cell_metadata %>%
    rename("barcode" = 1) %>%
    mutate(
        cell_type = str_replace_all(wg2_scpred_prediction, "_", " "),
        major_cell_type = case_when(
            wg2_scpred_prediction %in% c(
                "B_intermediate",
                "B_memory",
                "B_naive",
                "Plasmablast"
            ) ~ "B",
            wg2_scpred_prediction %in% c(
                "NK",
                "NK_CD56bright",
                "NK_Proliferating"
            ) ~ "NK",
            wg2_scpred_prediction %in% c(
                "CD8_Naive",
                "CD8_Proliferating",
                "CD8_TCM",
                "CD8_TEM"
            ) ~ "CD8 T",
            wg2_scpred_prediction %in% c(
                "CD4_CTL",
                "CD4_Naive",
                "CD4_Proliferating",
                "CD4_TCM",
                "CD4_TEM",
                "Treg"
            ) ~ "CD4 T",
            wg2_scpred_prediction %in% c(
                "dnT",
                "gdT",
                "ILC",
                "MAIT"
            ) ~ "Unconventional T",
            wg2_scpred_prediction %in% c(
                "pDC",
                "cDC1",
                "cDC2",
                "ASDC"
            ) ~ "Dendritic",
            wg2_scpred_prediction %in% c(
                "CD14_Mono",
                "CD16_Mono"
            ) ~ "Monocyte",
            wg2_scpred_prediction %in% c(
                "HSPC"
            ) ~ "Other"
        ) %>% fct_relevel(c("B", "NK", "CD8 T", "CD4 T", "Unconventional T", "Dendritic", "Monocyte", "Other"))
    )

# get the cell cycle phase data

# cell_cycle_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/240_libraries_cellcyle_phase.csv") %>%
#     rename("barcode" = 1)

# read in the UMAP coordinates
# umap_coords <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_harmony_umap_covs.csv") %>%
#     rename("barcode" = 1)

umap_coords <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_harmony_umap_coords_{data_type}.csv")) %>%
    rename("barcode" = 1)

plot_data <- umap_coords %>%
    left_join(cell_metadata) %>%
    # left_join(cell_cycle_metadata) %>%
    rename("UMAP 1" = UMAP1, "UMAP 2" = UMAP2)


# read in and combine cell cycle phase data

cell_cycle_csv_list <- list.files("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/per_library_phases/output", pattern = "_cell_cycle_scoring.csv", full.names = TRUE)

cell_cycle_data <- cell_cycle_csv_list %>%
    map(\(file) {
        read_csv(file = file)
    }) %>%
    list_rbind()

# ⚙️ Functions ----

square.ratio <- function(x) {
    range2 <- function(x) {
        sum(c(-1, 1) * range(x, na.rm = T))
    }
    return(coord_fixed(clip = "off", ratio = range2(x[, 1]) / range2(x[, 2])))
}

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

ggUMAPplot <- function(data,
                       group.by, # color by this metadata column
                       colorpal = "auto",
                       plot.title = NULL,
                       size = 0.1,
                       stroke = 0,
                       #    label = FALSE,
                       ...) {
    # data <- arrange(data, group.by)
    umap_ratio <- square.ratio(data[, c("UMAP 1", "UMAP 2")])

    new.plot <- ggplot(data) +
        geom_point(
            mapping = aes(
                x = .data[["UMAP 1"]],
                y = .data[["UMAP 2"]],
                colour = .data[[group.by]]
            ),
            size = size, stroke = stroke, ...
        ) +
        umap_ratio +
        theme_tenk10k() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            legend.text.align = 0,
            legend.key.size = unit(5, "mm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            ...
        ) +
        guides(color = guide_legend(override.aes = list(size = 4, shape = 15)))

    if (length(colorpal) > 1) {
        new.plot <- new.plot +
            scale_color_manual(values = colorpal, drop = T)
    }
    if (!is.null(plot.title)) {
        new.plot <- new.plot +
            labs(title = plot.title)
    } else {
        new.plot <- new.plot +
            labs(title = group.by)
    }

    # if (label) {
    #     Seurat::LabelClusters(new.plot, id = {{group.by}})
    # }

    return(new.plot)
}


# 📊 make the plots ----

# UMAP colored by cell type
plot_data <- plot_data %>%
    mutate(cell_type = factor(cell_type, levels = tenk_color_pal$cell_type))

cohort_plot <- plot_data %>%
    ggUMAPplot(group.by = "cohort")

cohort_plot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cohort.png"),
        width = 9, height = 5,
        dpi = 1300,
    )

cell_type_plot <- plot_data %>%
    ggUMAPplot(group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type))

cell_type_plot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cell_type.png"),
        width = 9, height = 5,
        dpi = 1300,
    )

# Umap colored by ell cycle phase

# plot_data %>%
#     ggUMAPplot(group.by = "phase") %>%
#     ggsave(
#         filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cell_cycle_phase.png",
#         width = 9, height = 5,
#         dpi = 1300
#     )

# Quality metric feature plots
# quality_metrics <- c("n_genes_by_counts", "total_counts", "pct_counts_mt", "S_score", "G2M_score")
quality_metrics <- c("n_genes_by_counts", "total_counts", "pct_counts_mt")

quality_metrics %>%
    purrr::walk(\(metric) ggFeaturePlot(
        data = plot_data,
        feature = metric,
        size = 0.1,
        save = TRUE,
        path = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/")
    ))


# Violin plots ----

# TODO: violin plot function

# ggViolin <- function(data) {

# }

# make violin plots for each cell type / qc metric

n_genes_by_counts_by_cell_types <- plot_data %>%
    pivot_longer(cols = {{ quality_metrics }}, names_to = "qc_metric", values_to = "qc_metric_value") %>%
    ggplot(aes(x = cell_type, y = qc_metric_value, fill = cell_type, color = cell_type)) +
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
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggVln_n_genes_by_counts_by_cell_types_300_libs.png"),
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

seq_lib_plot <- plot_data %>%
    ggUMAPplot(group.by = "sequencing_library")

seq_lib_plot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_seq_lib.png"),
        width = 10, height = 5,
        dpi = 1300,
    )
