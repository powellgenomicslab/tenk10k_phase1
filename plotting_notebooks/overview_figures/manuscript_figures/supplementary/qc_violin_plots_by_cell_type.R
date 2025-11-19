# conda activate r_python_dev

library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(data.table)


source("plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

# TODO: plot with nicely formatted names for quality metrics

# ⚙️ Parameters ----

# quality metrics to plot
quality_metrics <- c("n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo")
plot_titles <- c(
    "n_genes_by_counts" = "Number of genes detected",
    "total_counts" = "Total counts",
    "pct_counts_mt" = "Percentage mitochondrial counts",
    "pct_counts_ribo" = "Percentage ribosomal counts"
)

width <- 5
height <- 2
dpi <- 400

# 📗 Read in the data ----

# read in the cell-level metadata
cell_metadata <- get_latest_metadata()


# 📊 make the violin plots ----

plot_data <- cell_metadata %>%
    mutate(`Cell type` = factor(cell_type, levels = tenk_color_pal$cell_type))

n_genes_by_counts_by_cell_types <- plot_data %>%
    pivot_longer(cols = {{ quality_metrics }}, names_to = "qc_metric", values_to = "qc_metric_value") %>%
    mutate(qc_metric = factor(qc_metric, levels = quality_metrics)) %>%
    ggplot(aes(x = `Cell type`, y = qc_metric_value, fill = `Cell type`, color = `Cell type`)) +
    facet_wrap(~qc_metric, scales = "free_y", ncol = 1, labeller = as_labeller(plot_titles)) +
    geom_violin(width = 1) +
    scale_fill_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
    scale_color_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
    # guides(y = guide_axis(cap = "upper")) +
    # theme_classic() +
    theme_bw() + 
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        # panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none"
    )

n_genes_by_counts_by_cell_types %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/supplementary/figures/300_libraries_ggVln_qc_by_cell_types_no_cellbender.png"),
        width = width,
        height = height * length(quality_metrics),
        dpi = 300
    )
