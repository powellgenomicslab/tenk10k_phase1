library(data.table)
library(tidyverse)
library(patchwork)
library(ggsci)

# ⚙️ functions  ----

#' plot histogram of qc metric
qc_hist <- function(df, metric, ...) {
    df %>%
        ggplot(aes(x = !!sym(metric))) +
        geom_histogram(boundary = 0, fill = "skyblue", colour = "skyblue", bins = 200) +
        theme_classic() +
        coord_cartesian(xlim = c(min(df[[metric]], na.rm = TRUE), quantile(df[[metric]], 0.99, , na.rm = TRUE))) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), ...) %>%
        return()
}

qc_density <- function(df, metric) {
    df %>%
        ggplot(aes(x = !!sym(metric))) +
        geom_density() +
        theme_classic() +
        # coord_cartesian(xlim = c(min(df[[metric]], na.rm=TRUE),  quantile(df[[metric]], 0.98, , na.rm=TRUE))) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) %>%
        return()
}

qc_vln <- function(df, metric, group_col) {
    df %>%
        ggplot(aes(y = !!sym(metric), x = as.factor(!!sym(group_col)), fill = as.factor(!!sym(group_col)))) +
        geom_violin() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) %>%
        return()
}

# 📚 read in the data ----

tenk10k_metadata <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_cell_metadata_subset_pre_filtering.csv")
tenk10k_metadata <- tenk10k_metadata[, !duplicated(colnames(tenk10k_metadata)), with = FALSE] # oops. remove the duplicated columns

# compute the cell complexity score adata.obs['log10GenesPerUMI'] = np.log10(adata.obs['n_genes_by_counts']) / np.log10(adata.obs['total_counts'])
tenk10k_metadata[, log10_genes_per_umi := log10(n_genes_by_counts) / log10(total_counts), ]

# 📊 plot the data ----

metrics_to_plot <- c(
    "n_genes",
    "total_counts",
    "log10_genes_per_umi",
    "cellbender_background_fraction",
    "pct_counts_mt",
    "pct_counts_ribo",
    "pct_counts_hb"
)

hist_list <- metrics_to_plot %>%
    purrr::map(\(metric) qc_hist(df = tenk10k_metadata, metric = metric))

hist_list %>%
    patchwork::wrap_plots(ncol = 2) %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/qc/qc_histograms.png", height = 10, width = 8)

density_list <- metrics_to_plot %>%
    purrr::map(\(metric) qc_density(df = tenk10k_metadata, metric = metric))

density_list %>%
    patchwork::wrap_plots(ncol = 2) %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/qc/qc_density.png", height = 10, width = 8)

vln_list <- metrics_to_plot %>%
    purrr::map(\(metric) qc_vln(df = tenk10k_metadata, metric = metric, group_col = "leiden"))

vln_list %>%
    patchwork::wrap_plots(ncol = 2, guides = "collect") %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/qc/qc_vln_by_cluster_plots.png", height = 15, width = 15)

# zoomed in histograms ----

total_counts_hist <- qc_hist(df = tenk10k_metadata, metric = "total_counts")

total_counts_hist_zoomed <- qc_hist(df = tenk10k_metadata, metric = "total_counts") +
    geom_hline(yintercept = , colour = red, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 4000)) + scale_x_continuous(limits = c(0, 4000))

n_genes_hist <- qc_hist(df = tenk10k_metadata, metric = "n_genes")

n_genes_hist_zoomed <- qc_hist(df = tenk10k_metadata, metric = "n_genes") +
    geom_hline(yintercept = , colour = red, linetype = "dashed") +
    coord_cartesian(xlim = c(0, 2000)) + scale_x_continuous(limits = c(0, 2000))

hist_list_zoomed <- list(total_counts_hist, total_counts_hist_zoomed, n_genes_hist, n_genes_hist_zoomed)
hist_list_zoomed %>%
    patchwork::wrap_plots(ncol = 2, guides = "collect") %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/qc/qc_histograms_zoomed.png", height = 15, width = 15)

total_counts_hist_zoomed_split <- tenk10k_metadata %>%
    ggplot(aes(x = total_counts)) +
    geom_histogram(boundary = 0, fill = "skyblue", colour = "skyblue", bins = 200) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_continuous(limits = c(0, 4000)) +
    facet_wrap(~leiden, ncol = 2)

total_counts_hist_zoomed_split %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/figures/qc/qc_histograms_zoomed_split.png", height = 25, width = 15)


range(tenk10k_metadata$n_genes) # already a minimum threshold of 200
range(tenk10k_metadata$total_counts) # already a minimum threshold of 500

# examine some filtering thresholds ----

tenk10k_metadata <- tenk10k_metadata %>% mutate(
    # low_quality_cluster_or_platelets = if_else(leiden == 19 | wg2_scpred_prediction == "Platelet", "cluster19_or_platelet", "other"),
    low_quality_cluster = if_else(leiden == 19, "cluster19", "other"),
    below_1000_n_genes = if_else(n_genes < 1000, "remove", "keep"),
    below_800_n_genes = if_else(n_genes < 800, "remove", "keep"),
    below_700_n_genes = if_else(n_genes < 700, "remove", "keep"),
    below_600_n_genes = if_else(n_genes < 600, "remove", "keep"),
    below_500_n_genes = if_else(n_genes < 500, "remove", "keep"),
    below_1000_n_genes = if_else(n_genes < 500, "remove", "keep"),
    below_600_n_genes_or_below_800_total_counts = if_else(n_genes < 600 | total_counts < 800, "remove", "keep")
)

no_platelets <- tenk10k_metadata %>%
    filter(wg2_scpred_prediction != "Platelet")

table(no_platelets$below_600_n_genes_or_below_800_total_counts, no_platelets$low_quality_cluster)


# plot how many cells within and how many outside the problematic cluster get filtered with each threshold?

table(tenk10k_metadata$low_quality_cluster) # remove 21237 cells total
table(tenk10k_metadata$below_1000_n_genes, tenk10k_metadata$low_quality_cluster) # remove 88637 total cells
table(tenk10k_metadata$below_800_n_genes, tenk10k_metadata$low_quality_cluster) # remove 52852 total cells
table(tenk10k_metadata$below_700_n_genes, tenk10k_metadata$low_quality_cluster) # remove 39938 total cells - this one seems to remove most of the problematic cluster? with miniml collateral damage
table(tenk10k_metadata$below_600_n_genes, tenk10k_metadata$low_quality_cluster) # remove 28941 total cells
table(tenk10k_metadata$below_500_n_genes, tenk10k_metadata$low_quality_cluster)

tenk10k_metadata %>%
    ggplot() +
    geom_bar()
