# mm activate r-python-dev

library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(scales)
library(patchwork)
library(ggbeeswarm)

source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

# 📗 Read in the data ----

cell_metadata <- get_latest_metadata()

plot_data <- cell_metadata %>%
    select(cell_type, individual) %>%
    mutate(cell_type = factor(cell_type, levels = tenk_color_pal$cell_type))

# ⚙️ parameters ----

# set plot width, height, and dpi to be consistent across all the figures

width <- 5
height <- 3
dpi <- 400

# plotting functions go here

# 📊 make the plots ----

# get the plot order for cell types - applied to all plots.
# ordering by cell type frequency within each major cell type
cell_type_order <- cell_metadata %>%
    mutate(major_cell_type = fct_infreq(major_cell_type)) %>%
    group_by(major_cell_type) %>%
    mutate(cell_type = fct_infreq(cell_type)) %>%
    pull(cell_type) %>%
    levels()

barplot_ncells_per_celltype <- plot_data %>%
    ggplot(aes(x = fct_relevel(cell_type, cell_type_order), fill = cell_type)) +
    geom_bar() +
    labs(x = "Cell type", y = "Cells") +
    scale_fill_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = label_comma()) +
    theme_tenk10k() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")

barplot_ncells_per_celltype %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/barplot_ncells_per_celltype.png",
        width = width, height = height, dpi = dpi
    )

barplot_ncells_per_celltype_log10 <- plot_data %>%
    ggplot(aes(x = fct_relevel(cell_type, cell_type_order), fill = cell_type)) +
    geom_bar() +
    labs(x = "Cell type", y = bquote("Cells "(log[10]))) +
    scale_fill_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
    scale_y_continuous(trans = log10_trans(), expand = expansion(mult = c(0, 0.05)), labels = label_comma()) +
    theme_tenk10k() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")

barplot_ncells_per_celltype_log10 %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/barplot_ncells_per_celltype_log10.png",
        width = width, height = height, dpi = dpi
    )

# NOTE:
# currently this plot does not have data points for individuals if zero cells were detected for that cell type
# Do we still want to plot if the cell types in the individual was zero?

# For onek1k, this plot has the percentage of total cells for each individual on the y axis
vln_plot_ncells_per_individual_per_celltype <- plot_data %>%
    group_by(cell_type, individual) %>%
    count() %>%
    ungroup() %>%
    complete(cell_type, individual, fill = list("n" = 0)) %>%
    ggplot(aes(x = fct_relevel(cell_type, cell_type_order), y = n, fill = cell_type, colour = cell_type)) +
    geom_quasirandom(stroke = NA, size = 0.5) +
    # geom_boxplot(colour = "black", fill = "white") +
    labs(x = "Cell type", y = "Cells") +
    scale_fill_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
    scale_colour_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = label_comma()) +
    theme_tenk10k() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")
vln_plot_ncells_per_individual_per_celltype %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/vln_plot_ncells_per_individual_per_celltype.png",
        width = width, height = height, dpi = dpi
    )

# vln_plot_ncells_per_individual_per_celltype_pct <- plot_data %>%
#     group_by(individual) %>%
#     mutate(total_cells_individual = n()) %>%
#     group_by(cell_type, individual) %>%
#     mutate(total_cells_cell_type_individual = n()) %>%
#     reframe(percentage = total_cells_cell_type_individual / total_cells_individual) %>%
#     distinct() %>%
#     complete(cell_type, individual, fill = list(percentage = 0)) %>%
#     ggplot(aes(x = fct_relevel(cell_type, cell_type_order), y = percentage, fill = cell_type, colour = cell_type)) +
#     geom_quasirandom(stroke = NA, size = 0.5) +
#     # geom_violin(width=1.5) +
#     # geom_boxplot(colour = "black", fill = "white") +
#     labs(x = "Cell type", y = "Cells (%)")  +
#     scale_fill_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
#     scale_colour_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$cell_type)) +
#     scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = label_percent(accuracy = 1)) +
#     theme_tenk10k() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none")

# vln_plot_ncells_per_individual_per_celltype_pct %>%
#     ggsave(
#         filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/vln_plot_ncells_per_individual_per_celltype_pct.png",
#         width = width, height = height, dpi = dpi
#     )

ncells_per_individual_total_hist <- plot_data %>%
    group_by(individual) %>%
    count() %>%
    ggplot(aes(x = n)) +
    geom_histogram(binwidth = 250, boundary = 0, fill = "skyblue", colour = "skyblue") +
    labs(x = "Cells per individual", y = "Individuals") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = label_comma()) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)), labels = label_comma()) +
    theme_tenk10k() +
    theme(legend.position = "none")

ncells_per_individual_total_hist %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/ncells_per_individual_total_hist.png",
        width = width, height = height, dpi = dpi
    )

# combine the panels into a single image

combined_plots <- ncells_per_individual_total_hist +
    # barplot_ncells_per_celltype +
    barplot_ncells_per_celltype_log10 +
    vln_plot_ncells_per_individual_per_celltype +
    plot_layout(ncol = 1, axes = "collect")

combined_plots %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/figure_1/figures/fig1_bar_vln_combined.png",
        width = width + 1, height = height * 3, dpi = dpi
    )
