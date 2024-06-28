# conda activate r_python_dev

library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(Seurat)

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
        theme_bw() +
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

# create a nice color palette

# 📊 make the plots ----

# generate nested colors with material theme
red <- pal_material("red")(5)[2:5]
purple <- pal_material("purple")(4)[2:4]
pink <- pal_material("pink")(7)
deeppurple <- pal_material("deep-purple")(7)
indigo <- pal_material("indigo")(7)[]
blue <- pal_material("blue")(7)[]
lightblue <- pal_material("light-blue")(7)
teal <- pal_material("teal")(7)[]
lightgreen <- pal_material("light-green")(5)[2:5]
lime <- pal_material("lime")(7)[]
yellow <- pal_material("yellow")(7)[]
orange <- pal_material("orange")(7)[]
brown <- pal_material("brown")(4)[2:4]
grey <- pal_material("grey")(4)[2:4]

tenk_color_pal <- tribble(
    ~cell_type, ~color,
    # Lymphoid
    ## B cells
    "B_intermediate", red[1],
    "B_memory", red[4],
    "B_naive", orange[4],
    "Plasmablast", orange[7],
    ## NK cells
    "NK", pink[5],
    "NK_CD56bright", purple[1],
    "NK_Proliferating", purple[3],
    # CD8 T cells
    "CD8_Naive", deeppurple[3],
    "CD8_Proliferating", deeppurple[6],
    "CD8_TCM", indigo[4],
    "CD8_TEM", indigo[7],
    # CD4 T cells
    "CD4_CTL", blue[2],
    "CD4_Naive", blue[4],
    "CD4_Proliferating", blue[6],
    "CD4_TCM", lightblue[1],
    "CD4_TEM", lightblue[4],
    "Treg", lightblue[6],
    #
    "dnT", grey[3],
    "gdT", grey[6],
    "ILC", teal[3],
    "MAIT", teal[5],
    # myeloid
    ## DC
    "pDC", lime[3],
    "cDC1", lime[7],
    "cDC2", lightgreen[3],
    "ASDC", lightgreen[4],
    ## Monocyte
    "CD14_Mono", yellow[2],
    "CD16_Mono", yellow[6],
    #
    "HSPC", brown[3],
    #
    # "Platelet", brown[1],
    # "Eryth", brown[2],
    #
    # "Doublet", grey[2],
)

# UMAP colored by cell type
plot_data <- plot_data %>%
    mutate(wg2_scpred_prediction = factor(wg2_scpred_prediction, levels = tenk_color_pal$cell_type))

cell_types_plot <- plot_data %>%
    ggUMAPplot(group.by = "wg2_scpred_prediction", colorpal = setNames(tenk_color_pal$color, tenk_color_pal$cell_type))

cell_types_plot %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggumap_wg2_cell_types.png",
        width = 9, height = 5,
        dpi = 1300,
    )

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
