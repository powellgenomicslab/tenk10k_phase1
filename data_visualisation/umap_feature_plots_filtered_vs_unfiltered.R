# conda activate r_python_dev

library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(patchwork)
library(Seurat)
library(scales)

source("plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

data_type <- "filtered_min1000genes"

# 📗 Read in the data ----

# read in the cell-level metadata
# cell_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_cell_metadata_subset.csv") %>%
#     rename("barcode" = 1)

# cell_cycle_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/240_libraries_cellcyle_phase.csv") %>%
#     rename("barcode" = 1)

cell_metadata <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_{data_type}.csv")) %>%
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

# read in and combine cell cycle phase data

cell_cycle_csv_list <- list.files("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cell_cycle/per_library_phases/output", pattern = "_cell_cycle_scoring.csv", full.names = TRUE)

pattern <- ".*/(S\\d+[-\\w]*)_cell_cycle_scoring\\.csv"

# Extract sample names
# sample_names <- str_match(cell_cycle_csv_list, pattern)[,2]

cell_cycle_data <- cell_cycle_csv_list %>%
    map(\(file) {
        read_csv(file = file) %>%
            mutate(sequencing_library = str_match(file, pattern)[, 2]) %>%
            return()
    }) %>%
    list_rbind() %>%
    mutate(barcode = paste0(`...1` %>% str_remove(pattern = "-[0-9]+"), "_", sequencing_library)) %>%
    select(-`...1`, -sequencing_library)


plot_data <- umap_coords %>%
    left_join(cell_metadata) %>%
    left_join(cell_cycle_data, by = "barcode") %>%
    rename("UMAP 1" = UMAP1, "UMAP 2" = UMAP2) %>%
    mutate(cell_type = factor(cell_type, levels = tenk_color_pal$cell_type))

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
            # axis.ticks = element_blank(),
            # axis.text = element_blank(),
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


    return(new.plot)
}

# 📊 make the plots ----

# UMAP colored by cell type
# %>%
# mutate(problematic_cd14_mono = if_else(`UMAP 1` > -3 & `UMAP 1` < 0 & `UMAP 2` > 7 & `UMAP 2` < 11 & wg2_scpred_prediction == "CD14_Mono", TRUE, FALSE))


cohort_plot <- plot_data %>%
    ggUMAPplot(group.by = "cohort")

cohort_plot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cohort.png"),
        width = 9, height = 5,
        dpi = 1300,
    )

# filtered_problematic <- plot_data %>%
#     filter(`UMAP 1` > -3 & `UMAP 1` < 0 & `UMAP 2` > 7 & `UMAP 2` < 11) %>%
#     filter(wg2_scpred_prediction == "CD14_Mono")

# cohort_plot <- filtered_problematic %>%
#     ggUMAPplot(group.by = "cohort")

# cohort_plot %>%
#     ggsave(
#         filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cohort_problematic_cluster.png"),
#         width = 9, height = 5,
#         dpi = 1300,
#     )

cell_type_plot <- plot_data %>%
    ggUMAPplot(group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type))

cell_type_plot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cell_type.png"),
        width = 9, height = 5,
        dpi = 1300,
    )


cell_typist_plot <- plot_data %>%
    ggUMAPplot(group.by = "celltypist_majority_voting")

cell_typist_plot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_celltypist_majority_voting.png"),
        width = 9, height = 5,
        dpi = 1300,
    )

# individuals_plot <- filtered_problematic %>%
#     ggUMAPplot(group.by = "individual")

# individuals_plot %>%
#     ggsave(
#         filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_individual_filtered_problematic.png"),
#         width = 9, height = 5,
#         dpi = 1300,
#     )

# problematic_libs <- paste0("S0", 167:210)

# plot_data <- plot_data %>%
#     mutate(problematic_library = if_else(sequencing_library %in% problematic_libs, "Outlier library", "Not outlier library"))

# only_cd14_mono_problematic_libraries <- plot_data %>%
#     filter(wg2_scpred_prediction == "CD14_Mono") %>%
#     filter(problematic_library == "Outlier library")

# only_cd14_mono_problematic_libraries %>%
#     group_by(individual, problematic_cd14_mono) %>%
#     count() %>%
#     pivot_wider(names_from = problematic_cd14_mono, values_from = n) %>%
#     filter(`TRUE` > `FALSE`) %>%
#     pull(individual)

# table(only_cd14_mono_new_pools$problematic_cd14_mono, only_cd14_mono_new_pools$individual)

# problematic_by_indiv_bar <- only_cd14_mono_problematic_libraries %>%
#     ggplot(aes(x = individual, groups = problematic_cd14_mono, fill = problematic_cd14_mono)) +
#     facet_wrap(~sequencing_library, nrow = 1, scales = "free_x") +
#     geom_bar() +
#     labs(title = "CD14 Mono") +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# problematic_by_indiv_bar %>%
#     ggsave(
#         filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_filtered_bar_seq_lib_problematic_cluster.png",
#         width = 25, height = 5
#     )

# Umap colored by ell cycle phase
# plot_data %>%
#     ggUMAPplot(group.by = "phase") %>%
#     ggsave(
#         filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cell_cycle_phase.png"),
#         width = 9, height = 5,
#         dpi = 1300
#     )

# Quality metric feature plots ----
# quality_metrics <- c("n_genes_by_counts", "total_counts", "pct_counts_mt", "S_score", "G2M_score")
quality_metrics <- c("n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "cellbender_background_fraction")

quality_metrics %>%
    purrr::walk(\(metric) ggFeaturePlot(
        data = plot_data,
        feature = metric,
        size = 0.1,
        save = TRUE,
        path = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/")
    ))

cell_cycle_metrics <- c("S_score", "G2M_score")

summary(plot_data$S_score)
summary(plot_data$G2M_score)

plot_list <- cell_cycle_metrics %>%
    purrr::map(\(metric) ggFeaturePlot(
        data = plot_data,
        feature = metric,
        size = 0.1,
        save = FALSE
    ))

plot_list %>%
    patchwork::wrap_plots(ncol = 2) %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_cell_cycle_phase_scores.png"), width = 15, height = 5)

# Violin plots ----

# make violin plots for each cell type / qc mWetric

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
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggVln_qc_by_cell_typess.png"),
        width = 14,
        height = 14,
        dpi = 300
    )

# TODO: make violin plots for each sequencing batch

# export ncells per sample

# cell_metadata %>%
#     group_by(cpg_id, cohort) %>%
#     summarise(total_barcodes = n()) %>%
#     write_csv("240_lib_all_individuals_ncells.csv")

# seq_lib_plot <- filtered_problematic %>%
#     ggUMAPplot(group.by = "sequencing_library")

# seq_lib_plot %>%
#     ggsave(
#         filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_seq_lib_problematic_cluster.png"),
#         width = 15, height = 5,
#         dpi = 1300,
#     )

# problematic_bar <- filtered_problematic %>%
#     filter(cell_type %in% c("CD14 Mono", "pDC")) %>%
#     ggplot(aes(x = sequencing_library)) +
#     geom_bar() +
#     facet_wrap(~cell_type, ncol = 1) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# problematic_bar %>%
#     ggsave(
#         filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_bar_seq_lib_problematic_cluster.png"),
#         width = 15, height = 10,
#         dpi = 1300,
#     )

# problematic_libs <- paste0("S0", 167:210)

# plot_data <- plot_data %>%
#     mutate(problematic_library = if_else(sequencing_library %in% problematic_libs, "Outlier library", "Not outlier library"))

# problematic_only_umap <- plot_data %>%
#     filter(problematic_library == "Outlier library") %>%
#     ggUMAPplot(group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type))

# problematic_only_umap %>%
#     ggsave(
#         filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cell_type_problematic_libraries_only.png"),
#         width = 10, height = 5,
#         dpi = 1300,
#     )

# no_problematic_umap <- plot_data %>%
#     filter(problematic_library != "Outlier library") %>%
#     ggUMAPplot(group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type))

# no_problematic_umap %>%
#     ggsave(
#         filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cell_type_no_problematic.png"),
#         width = 10, height = 5,
#         dpi = 1300,
#     )

# plot_data %>%
#     filter(sequencing_library == "S0041") %>%
#     pull(cpg_id) %>%
#     unique()



cell_cycle_umap <- ggplot(plot_data) +
    geom_point(
        mapping = aes(
            x = `UMAP 1`,
            y = `UMAP 2`,
            colour = phase,
        ),
        size = 0.1,
        stroke = 0,
        alpha = 0.2
    ) +
    theme_tenk10k() +
    theme(
        aspect.ratio = 1,
        legend.position = "right",
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.key.size = unit(5, "mm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()
    ) +
    guides(color = guide_legend(override.aes = list(size = 4, shape = 15))) +
    scale_color_manual(values = c("G1" = "lightgrey", "G2M" = "yellow", "S" = "red"), drop = T) +
    facet_wrap(~phase)

# Umap colored by ell cycle phase
cell_cycle_umap %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/300_libraries_{data_type}_ggumap_cell_cycle_phase_modified.png"),
        width = 9, height = 5,
        dpi = 1300
    )

# only_cd14_mono_problematic_libraries %>%
#     group_by(ct_id, problematic_cd14_mono) %>%
#     count() %>%
#     pivot_wider(names_from = problematic_cd14_mono, values_from = n) %>%
#     mutate(fraction = `TRUE` / (`FALSE` + `TRUE`)) %>%
#     filter(fraction > 0.8) %>%
#     arrange(desc(fraction)) %>%
#     select(ct_id) %>%
#     print(n = 100)


# highfraction_ind <- c(
#     "CPG508077",
#     "CPG499004",
#     "CPG498931",
#     "CPG500686",
#     "CPG498949",
#     "CPG501205",
#     "CPG500553",
#     "CPG500777",
#     "CPG499046",
#     "CPG499038",
#     "CPG501619",
#     "CPG498972",
#     "CPG498956",
#     "CPG499210",
#     "CPG508051",
#     "CPG498998",
#     "CPG499012",
#     "CPG500579",
#     "CPG499780",
#     "CPG501007"
# )

# only_cd14_mono_problematic_libraries %>%
#     filter(individual %in% highfraction_ind) %>%
#     select(ct_id) %>%
#     distinct()
