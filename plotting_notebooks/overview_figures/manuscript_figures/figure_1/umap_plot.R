library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(data.table)

source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

# 📗 Read in the data ----

# read in the cell-level metadata

cell_metadata <- get_latest_metadata()

# read in the UMAP coordinates

umap_coords <- get_latest_umap()

plot_data <- umap_coords %>%
    left_join(cell_metadata) %>%
    rename("UMAP 1" = UMAP1, "UMAP 2" = UMAP2) %>%
    mutate(cell_type = factor(cell_type, levels = tenk_color_pal$cell_type))

# ⚙️ Functions ----

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

# UMAP colored by cell type for figure 1

cell_types_plot <- plot_data %>%
    ggUMAPplot(group.by = "cell_type", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type), plot.title = NA)

cell_types_plot %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/umap_wg2_cell_types_dpi_1300.png",
        width = 9, height = 5,
        dpi = 1300,
    )

cell_types_plot %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/umap_wg2_cell_types_dpi_500.png",
        width = 9, height = 5,
        dpi = 500,
    )
