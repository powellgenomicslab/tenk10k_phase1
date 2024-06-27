library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(patchwork)
library(tidyverse)

# read in the cell-level metadata
cell_metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_archive/240_libraries_celltypes_by_cohort_indiv.csv") %>%
    rename("barcode" = 1)

cell_metadata %>%
    colnames()

plot_data <- cell_metadata %>%
    group_by(wg2_scpred_prediction, individual) %>%
    # calculate which samples have cells above each threshold
    mutate(
        three_or_more_cells = n() >= 3,
        five_or_more_cells = n() >= 5,
        ten_or_more_cells = n() >= 10,
        twenty_or_more_cells = n() >= 20,
        fifty_or_more_cells = n() >= 50,
        total_cells = n()
    ) %>%
    ungroup() %>%
    select(individual, wg2_scpred_prediction, three_or_more_cells, five_or_more_cells, ten_or_more_cells, twenty_or_more_cells, fifty_or_more_cells, total_cells) %>%
    distinct() %>%
    arrange(total_cells) %>%
    filter(three_or_more_cells == TRUE)

bar_plot_passing_thresholds <- plot_data %>%
    group_by(wg2_scpred_prediction) %>%
    mutate(total_individuals_with_3plus_cells = sum(three_or_more_cells)) %>%
    ungroup() %>%
    select(wg2_scpred_prediction, total_individuals_with_3plus_cells) %>%
    distinct() %>%
    ggplot(aes(x = reorder(wg2_scpred_prediction, -total_individuals_with_3plus_cells), y = total_individuals_with_3plus_cells)) +
    geom_bar(fill = "skyblue", stat = "identity") +
    theme_bw() +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "red") +
    xlab("Cell type") +
    ylab("Number of individuals with 3 or more cells") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

bar_plot_passing_thresholds %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggBarplot_cell_types_passing_thresholds.png",
        width = 16,
        height = 7,
        dpi = 300
    )


bar_plot_total_cells <- plot_data %>%
    group_by(wg2_scpred_prediction) %>%
    summarize(total_n_cells = sum(total_cells)) %>%
    ggplot(aes(x = reorder(wg2_scpred_prediction, -total_n_cells), y = total_n_cells)) +
    geom_bar(fill = "skyblue", stat = "identity") +
    theme_bw() +
    xlab("Cell type") +
    ylab("Total number of cells") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

bar_plot_total_cells %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggBarplot_total_cells.png",
        width = 16,
        height = 7,
        dpi = 300
    )

combined <- bar_plot_passing_thresholds + bar_plot_total_cells + plot_layout(ncol = 1, axes = "collect")

combined %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggBarplot_over_3_cells_combined_plots.png",
        width = 10,
        height = 8,
        dpi = 300
    )

cell_threshold_data <- cell_metadata %>%
    group_by(wg2_scpred_prediction, individual) %>%
    # calculate which samples have cells above each threshold
    mutate(
        three_or_more_cells = n() >= 3,
        five_or_more_cells = n() >= 5,
        ten_or_more_cells = n() >= 10,
        twenty_or_more_cells = n() >= 20,
        fifty_or_more_cells = n() >= 50,
        total_cells = n()
    ) %>%
    ungroup() %>%
    select(individual, wg2_scpred_prediction, three_or_more_cells, five_or_more_cells, ten_or_more_cells, twenty_or_more_cells, fifty_or_more_cells) %>%
    distinct() %>%
    group_by(wg2_scpred_prediction) %>%
    summarise(
        three_or_more_cells = sum(three_or_more_cells),
        five_or_more_cells = sum(five_or_more_cells),
        ten_or_more_cells = sum(ten_or_more_cells),
        twenty_or_more_cells = sum(twenty_or_more_cells),
        fifty_or_more_cells = sum(fifty_or_more_cells)
    ) %>%
    ungroup() %>%
    pivot_longer(!wg2_scpred_prediction, names_to = "cell_number_threshold", values_to = "total_individuals") %>%
    mutate(cell_number_threshold = factor(cell_number_threshold, levels = c(
        "three_or_more_cells",
        "five_or_more_cells",
        "ten_or_more_cells",
        "twenty_or_more_cells",
        "fifty_or_more_cells"
    )))

bar_plot_passing_thresholds_multi <- cell_threshold_data %>%
    ggplot(aes(x = reorder(wg2_scpred_prediction, -total_individuals), y = total_individuals)) +
    facet_wrap(~cell_number_threshold, ncol = 1) +
    geom_bar(fill = "skyblue", stat = "identity") +
    theme_bw() +
    geom_hline(yintercept = 100, linetype = "dashed", colour = "red") +
    xlab("Cell type") +
    ylab("Number of individuals") +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

bar_plot_passing_thresholds_multi %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/ggBarplot_cell_types_passing_thresholds_multi_thresholds.png",
        width = 14,
        height = 16,
        dpi = 300
    )
