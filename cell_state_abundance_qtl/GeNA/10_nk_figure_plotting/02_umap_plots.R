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
    q1 <- data %>%
        pull({{ feature }}) %>%
        quantile(0.01, na.rm = TRUE)
    q99 <- data %>%
        pull({{ feature }}) %>%
        quantile(0.99, na.rm = TRUE)

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
        theme_bw() +
        theme(
            aspect.ratio = 1,
            legend.position = "right",
            legend.key.size = unit(5, "mm"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            plot.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            ...
        ) +
        # scale_color_gradientn(colours = c("lightgrey", "red"))
        scale_color_viridis(
            option = "magma",
            direction = -1,
            limits = c(q1, q99),
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

