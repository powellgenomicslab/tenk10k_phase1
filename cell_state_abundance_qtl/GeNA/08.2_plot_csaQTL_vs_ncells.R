library(tidyverse)
library(data.table)
library(glue)
library(scales)
library(ggsci)

source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt")
celltypes <- celltypes[!celltypes %in% c("ALL", "Other")] # remove Other as there were no significant csaqtls

# GeNA GWAS summary statistics
read_summstats <- function(celltype, analysis_name, resolution) {
    file_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats_lead_snps_MAF_0.05.tsv")
    sumstats <- fread(file_path, select = 1:10) # read in results without the beta's
    sumstats[, `:=`(
        P = fifelse(as.numeric(P) == 0, .Machine$double.xmin * .Machine$double.eps, as.numeric(P)),
        celltype = celltype
    )]
    return(sumstats)
}

# combine summary statistics for all cell types
sumstats_all_ct <- celltypes %>%
    map(\(celltype) read_summstats(celltype = celltype, analysis_name = analysis_name, resolution = resolution)) %>%
    rbindlist(fill = TRUE)
sumstats_all_ct[, celltype := factor(str_replace(celltype, pattern = "_", replacement = " "), levels = unique(tenk_color_pal$major_cell_type)), ]

sumstats_all_ct_sig_summary <- sumstats_all_ct[P < 5e-8, .N, by = celltype]

# save combined file for lead snps
sumstats_all_ct %>% fwrite("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/major_cell_types/summstats_all_lead_csaqtl.csv")

# number of unique csaQTL loci (number of unique lead SNPs)
sumstats_all_ct$ID %>%
    n_distinct()

metadata <- get_latest_metadata()

ncells_nsnps_per_ct <- metadata %>%
    group_by(major_cell_type) %>%
    count() %>%
    rename("celltype" = major_cell_type, "Number of cells in cell type" = `n`) %>%
    left_join(sumstats_all_ct_sig_summary, by = "celltype") %>%
    rename("Number of significant csaQTLs" = N) %>%
    mutate(celltype = factor(celltype, levels = unique(tenk_color_pal$major_cell_type))) %>%
    mutate(
        `Number of significant csaQTLs` = replace_na(`Number of significant csaQTLs`, replace = 0)
    )


ncells_by_nsnps_plot <- ncells_nsnps_per_ct %>%
    ggplot(aes(y = `Number of significant csaQTLs`, x = `Number of cells in cell type`, colour = celltype)) +
    geom_point(size = 3) +
    theme_classic() +
    scale_x_continuous(labels = label_comma()) +
    scale_y_continuous(labels = label_comma()) +
    scale_colour_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    theme(aspect.ratio = 1)

ncells_by_nsnps_plot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/combined_ncells_by_n_lead_snps_plot_{analysis_name}.png"),
        width = 6, height = 5
    )

ncells_nsnps_per_ct %>%
    write_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/{resolution}/ncells_n_lead_snps_per_ct_{analysis_name}.csv"))

lead_csaqtl_bar <- ncells_nsnps_per_ct %>%
    ggplot(aes(x =  fct_reorder(celltype, `Number of significant csaQTLs`, .desc = TRUE, .na_rm = FALSE), y = `Number of significant csaQTLs`, colour = celltype, fill = celltype)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    # scale_x_continuous(labels = label_comma()) +
    scale_y_continuous(labels = label_comma()) +
    scale_colour_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    scale_fill_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    theme(aspect.ratio = 0.5) +
    labs(x = "Cell type") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

lead_csaqtl_bar %>%
    ggsave(filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/n_lead_snps_per_ct_bar_{analysis_name}.png"), width = 8, height = 4)
