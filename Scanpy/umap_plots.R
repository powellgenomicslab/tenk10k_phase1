library(tidyverse)
library(ggsci)

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
                       size = 0.001,
                       stroke = 0, ...) {
    data <- arrange(data, group.by)
    umap.ratio <- square.ratio(data[, c("UMAP1", "UMAP2")])

    new.plot <- ggplot(data) +
        geom_point(
            mapping = aes(x = .data[["UMAP1"]], y = .data[["UMAP2"]], colour = .data[[group.by]]),
            size = size, stroke = stroke, ...
        ) +
        umap.ratio +
        theme_bw() +
        theme(
            legend.position = "right",
            legend.title = element_blank(),
            legend.text.align = 0,
            legend.key.size = unit(5, "mm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
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


# read in the cell-level metadata
cell_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_cell_metadata_subset.csv") %>%
    rename("barcode" = 1)

# read in the UMAP coordinates
umap_coords <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_harmony_umap_covs.csv") %>%
    rename("barcode" = 1)


plot_data <- umap_coords %>%
    left_join(cell_metadata)

head(cell_metadata)
dim(cell_metadata)
colnames(cell_metadata)

ggUMAPplot(data = plot_data, group.by = "wg2_scpred_prediction") %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggumap_wg2_cell_types.png",
        width = 9, height = 5,
        dpi = 1000
    )
