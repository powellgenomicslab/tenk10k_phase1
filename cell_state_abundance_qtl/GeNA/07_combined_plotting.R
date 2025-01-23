library(tidyverse)
library(qqman)
library(data.table)
library(glue)
library(scattermore)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

analysis_name <- "with_expr_pc_covars"

# 📚 read in the data ----

celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/unique_cell_types_wg2_scpred.txt")

read_summstats <- function(celltype, analysis_name) {
    sumstats <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{celltype}/{analysis_name}/GeNA_sumstats.txt")) %>%
        .[, P := as.numeric(P), ] %>%
        .[, P := ifelse(P == 0, .Machine$double.xmin * .Machine$double.eps, P)] %>% # replace P of zero with minimum floating point value
        .[, celltype := celltype] %>%
        .[, grep("BETA_NAMPC", names(.)) := NULL, ] # remove betas 
    return(sumstats)
}
# combine summary statistics for all cell types
sumstats_all_ct <- celltypes[!celltypes == "NK"] %>%
    map(\(celltype) read_summstats(celltype = celltype, analysis_name = analysis_name)) %>%
    rbindlist(fill = TRUE)

sumstats_all_ct[,celltype := factor(celltype, levels = tenk_color_pal$wg2_scpred_prediction),]

# remove the betas
# sumstats_all_ct[, grep("BETA_NAMPC", names(sumstats_all_ct)) := NULL, ]

setorder(sumstats_all_ct, P)

sumstats_all_ct[
    ,
    observed := -log10(P),
]

# calculate -log10 expected p value for each cell type 
sumstats_all_ct[
    ,
    expected := -log10(ppoints(.N)),
    by = celltype
]

log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

combined_qqplot <- sumstats_all_ct %>%
    ggplot(aes(x = expected, y = observed, colour = celltype)) +
    labs(
        x = log10Pe,
        y = log10Po
    ) +
    # geom_point() + # can try scattermost for faster?
    geom_scattermore() + 
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    theme_classic() +
    scale_colour_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$wg2_scpred_prediction))

# add colour palette

combined_qqplot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/combined_plots/combined_qq_{analysis_name}.png"),
        width = 10, height = 7
    )


#-----------------------------------------------
# Number of cells in cell type vs number of significant snps?
#-----------------------------------------------


# https://www.nature.com/articles/s41588-024-01909-1#:~:text=Therefore%2C%20we%20consider%20loci%20with%20GeNA%20P%E2%80%89%3C%E2%80%895%E2%80%89%C3%97%E2%80%8910%E2%88%928%20associations%20genome%2Dwide%20significant.
# GeNA manuscript: "we consider loci with GeNA P < 5 × 10−8 associations genome-wide significant."

sumstats_all_ct_sig_summary <- sumstats_all_ct[P < 5e-8, .N, by = celltype]

metadata <- get_latest_metadata()

ncells_nsnps_per_ct <- metadata %>%
    group_by(wg2_scpred_prediction) %>%
    count() %>%
    rename("celltype" = wg2_scpred_prediction, "Number of cells in cell type" = `n`) %>%
    left_join(sumstats_all_ct_sig_summary, by = "celltype") %>%
    rename("Number of significant csaQTLs" = N)

ncells_by_nsnps_plot <- ncells_nsnps_per_ct %>%
    ggplot(aes(y = `Number of significant csaQTLs`, x = `Number of cells in cell type`, colour = celltype)) +
    geom_point() +
    theme_classic() +
    scale_colour_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$wg2_scpred_prediction))

ncells_by_nsnps_plot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/combined_plots/combined_ncells_by_nsnps_plot_{analysis_name}.png"),
        width = 10, height = 7
    )

ncells_nsnps_per_ct %>% write_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/ncells_nsnps_per_ct.csv")

# subset to just one SNP per LD-block? i.e. just the lead SNPs?
