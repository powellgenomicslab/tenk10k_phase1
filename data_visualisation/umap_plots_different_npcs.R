library(tidyverse)
library(glue)
library(patchwork)

source("plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

# ⚙️ plotting functions ----

square.ratio <- function(x) {
    range2 <- function(x) {
        sum(c(-1, 1) * range(x, na.rm = T))
    }
    return(coord_fixed(clip = "off", ratio = range2(x[, 1]) / range2(x[, 2])))
}

ggUMAPplot <- function(data,
                       group.by, # color by this metadata column
                       colorpal = "auto",
                       plot.title = NULL,
                       umap_col_names = c("UMAP 1", "UMAP 2"),
                       size = 0.1,
                       stroke = 0,
                       #    label = FALSE,
                       ...) {
    # data <- arrange(data, group.by)
    umap_ratio <- square.ratio(data[, umap_col_names])

    new.plot <- ggplot(data) +
        geom_point(
            mapping = aes(
                x = .data[[umap_col_names[1]]],
                y = .data[[umap_col_names[2]]],
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

# UMAP with different number of principal components ----

# get the predicted cell type data
cell_metadata <- get_latest_metadata() %>%
    select(barcode, cell_type, cohort) %>%
    mutate(cell_type = factor(cell_type, levels = tenk_color_pal$cell_type))

# get the current UMAP (50 PCs)
umap_coords_50_pcs <- get_latest_umap() %>%
    rename(
        "UMAP_1_50_PCs" = UMAP1,
        "UMAP_2_50_PCs" = UMAP2,
    )

umap_coords_list <- list()
umap_coords_list[["UMAP_50_PCs"]] <- umap_coords_50_pcs

pc_list <- seq(10, 45, 5)
pc_list <- c(30)

# get the UMAP generated with different number of PCs
# TODO: add in 45 PC's when it's ready
for (n_pcs in pc_list) {
    umap_1_col_name <- glue("UMAP_1_{n_pcs}_PCs")
    umap_2_col_name <- glue("UMAP_2_{n_pcs}_PCs")

    umap_coords <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/umap_coords_{n_pcs}_pcs.csv")) %>%
        dplyr::rename(
            !!umap_1_col_name := UMAP1,
            !!umap_2_col_name := UMAP2,
            barcode = 1
        )

    item_name <- glue("UMAP_{n_pcs}_PCs")
    umap_coords_list[[item_name]] <- umap_coords
}

# join together each UMAP with the cell types data
plot_data <- umap_coords_list %>%
    reduce(left_join, by = "barcode") %>%
    left_join(cell_metadata, by = "barcode")

# individual plot for 30 PC's

cell_types_plot <- plot_data %>%
    ggUMAPplot(umap_col_names = c("UMAP_1_30_PCs", "UMAP_2_30_PCs"), group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type), plot.title = NA)

cell_types_plot %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/umap_wg2_cell_types_30_pcs_new_cols.png",
        width = 9, height = 5,
        dpi = 500,
    )

plot_list <- list()
for (n_pcs in c(pc_list, 50)) {
    cell_types_plot <- plot_data %>%
        ggUMAPplot(
            umap_col_names = c(glue("UMAP_1_{n_pcs}_PCs"), glue("UMAP_2_{n_pcs}_PCs")),
            group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type), plot.title = glue("{n_pcs} principal components")
        )

    cell_types_plot %>%
        ggsave(
            filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/umap_wg2_cell_types_{n_pcs}_pcs.png"),
            width = 9, height = 5,
            dpi = 500,
        )

    plot_list[[as.character(n_pcs)]] <- cell_types_plot
}

combined_umaps <- wrap_plots(plot_list) +
    plot_layout(guides = "collect", ncol = 3)

combined_umaps %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/umap_wg2_cell_types_different_pcs_combined.png"),
        width = 18, height = 12,
        dpi = 500,
    )

# UMAP with different number of HVG ----

hvg_list <- seq(500, 3000, 500)

umap_coords_list <- list()

for (n_hvgs in hvg_list) {
    umap_1_col_name <- glue("UMAP_1_{n_hvgs}_HVGs")
    umap_2_col_name <- glue("UMAP_2_{n_hvgs}_HVGs")

    umap_coords <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/umap_coords_{n_hvgs}_hvg.csv")) %>%
        dplyr::rename(
            !!umap_1_col_name := UMAP1,
            !!umap_2_col_name := UMAP2,
            barcode = 1
        )

    item_name <- glue("UMAP_{n_hvgs}_HVGs")
    umap_coords_list[[item_name]] <- umap_coords
}

plot_data <- umap_coords_list %>%
    reduce(left_join, by = "barcode") %>%
    left_join(cell_metadata, by = "barcode")

plot_list <- list()
for (n_hvgs in hvg_list) {
    cell_types_plot <- plot_data %>%
        ggUMAPplot(
            umap_col_names = c(glue("UMAP_1_{n_hvgs}_HVGs"), glue("UMAP_2_{n_hvgs}_HVGs")),
            group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type), plot.title = glue("{n_hvgs} highly variable genes")
        )

    cell_types_plot %>%
        ggsave(
            filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/umap_wg2_cell_types_{n_hvgs}_hvgs.png"),
            width = 9, height = 5,
            dpi = 500,
        )

    plot_list[[as.character(n_hvgs)]] <- cell_types_plot
}

combined_umaps <- wrap_plots(plot_list) +
    plot_layout(guides = "collect", ncol = 3)

combined_umaps %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/umap_wg2_cell_types_different_hvgs_combined.png"),
        width = 18, height = 12,
        dpi = 500,
    )

# UMAP for specific cohort and combined cohort - note these are the only ones that have been run with the outlier cell type composition samples removed ----

cohort_list <- c("BioHEART", "TOB")
umap_coords_list <- list()

# get plot data for each cohort 
for (cohort in cohort_list) {
    umap_1_col_name <- glue("UMAP_1_3000_HVGs_30_PCs_{cohort}")
    umap_2_col_name <- glue("UMAP_2_3000_HVGs_30_PCs_{cohort}")

    umap_coords <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/umap_coords_30_pcs_3000_hvg_{cohort}_only_outliers_removed.csv")) %>%
        dplyr::rename(
            !!umap_1_col_name := UMAP1,
            !!umap_2_col_name := UMAP2,
            barcode = 1
        )

    item_name <- glue("UMAP_{cohort}_cohort")
    umap_coords_list[[item_name]] <- umap_coords
}
# manually add in the combined cohort (file name doesnt work with for-loop above). This version has the outlier samples removed 
umap_coords_list[["UMAP_combined_cohort"]] <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/umap_coords_30_pcs_3000_hvg_outliers_removed.csv")) %>%
    dplyr::rename(
        "UMAP_1_3000_HVGs_30_PCs_combined" := UMAP1,
        "UMAP_2_3000_HVGs_30_PCs_combined" := UMAP2,
        barcode = 1
    )

plot_data <- umap_coords_list %>%
    reduce(full_join, by = "barcode") %>%
    left_join(cell_metadata, by = "barcode")

plot_list <- list()
for (cohort in c(cohort_list, "combined")) {
    cell_types_plot <- plot_data %>%
        ggUMAPplot(
            umap_col_names = c(glue("UMAP_1_3000_HVGs_30_PCs_{cohort}"), glue("UMAP_2_3000_HVGs_30_PCs_{cohort}")),
            group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type), plot.title = glue("{cohort} cohort")
        )

    cell_types_plot %>%
        ggsave(
            filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/umap_wg2_cell_types_30_pcs_3000_hvg_{cohort}_outliers_removed.png"),
            width = 9, height = 5,
            dpi = 500,
        )

    plot_list[[as.character(cohort)]] <- cell_types_plot
}

plot_list[["umap_combined_cohort_colors"]] <-  plot_data %>%
        ggUMAPplot(
            umap_col_names = c(glue("UMAP_1_3000_HVGs_30_PCs_combined"), glue("UMAP_2_3000_HVGs_30_PCs_combined")),
            group.by = "cohort", colorpal = c("#E41A1C", "#377EB8"), plot.title = glue("{cohort} cohort")
        )


# make the combined plot 
combined_umaps <- wrap_plots(plot_list) +
    plot_layout(guides = "collect", ncol = 2)

combined_umaps %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/umap_wg2_cell_types_separate_cohorts_outliers_removed_combined.png"),
        width = 18, height = 12,
        dpi = 500,
    )
