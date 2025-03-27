library(tidyverse)
library(viridis)
library(scales)
library(glue)
library(scattermore)
library(ggbeeswarm)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")


# ⚙️ Functions ----

ggUMAPplot <- function(
    data,
    group.by, # color by this metadata column
    colorpal = "auto",
    plot.title = NULL,
    size = 0.1,
    stroke = 0,
    ptsize = 2.5,
    px = 350,
    #    label = FALSE,
    ...) {
    new.plot <- ggplot(data) +
        geom_scattermore(
            mapping = aes(x = .data[["UMAP 1"]], y = .data[["UMAP 2"]], colour = .data[[group.by]]),
            pointsize = ptsize, pixels = c(px * 5, px * 5), ...
        ) +
        theme_tenk10k() +
        theme(
            aspect.ratio = 1,
            legend.position = "right",
            legend.title = element_blank(),
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

ggFeaturePlot <- function(
    data,
    feature,
    save = FALSE,
    path = NULL,
    size = 0.1,
    stroke = 0,
    ptsize = 2.5,
    px = 350,
    ...) {
    data <- arrange(data, feature) # plot the lowest expressing cells first and highest expression last

    # scale colours between 1st and 99th quantile to reduce impact of outliers on the colour scaling
    qmin <- data %>%
        pull({{ feature }}) %>%
        quantile(0.05, na.rm = TRUE)
    qmax <- data %>%
        pull({{ feature }}) %>%
        quantile(0.95, na.rm = TRUE)

    # umap_ratio <- square.ratio(data[, c("UMAP 1", "UMAP 2")])
    feature_plot <- ggplot(data) +
        # geom_point(
        #     mapping = aes(x = .data[["UMAP 1"]], y = .data[["UMAP 2"]], colour = .data[[feature]]),
        #     size = size, stroke = stroke, ...
        # ) +
        geom_scattermore(
            mapping = aes(x = .data[["UMAP 1"]], y = .data[["UMAP 2"]], colour = .data[[feature]]),
            pointsize = ptsize, pixels = c(px * 5, px * 5), ...
        ) +
        labs(title = feature) +
        theme_tenk10k() +
        theme(
            aspect.ratio = 1,
            legend.position = "right",
            legend.title = element_blank(),
            legend.key.size = unit(5, "mm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            ...
        ) +
        # scale_color_gradientn(colours = c("lightgrey", "red"))
        scale_color_viridis(
            option = "magma",
            direction = 1,
            begin = 0,
            limits = c(qmin, qmax),
            oob = scales::squish
        )

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
