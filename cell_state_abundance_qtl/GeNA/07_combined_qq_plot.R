library(tidyverse)
library(data.table)
library(glue)
library(scattermore)
library(scales)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

# 📚 read in the data ----

# celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/unique_cell_types_wg2_scpred.txt")
celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt")
celltypes <- celltypes[celltypes != "ALL"]

# Get the minor allele frequencies
read_afreq <- function(afreq_path) {
    read_tsv(
        afreq_path,
        comment = "#",
        col_names = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"),
        col_types = cols(
            CHROM = col_character(),
            ID = col_character(),
            REF = col_character(),
            ALT = col_character(),
            ALT_FREQS = col_double(),
            OBS_CT = col_integer(),
        )
    )
}
afreq_path <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15.afreq"
afreq <- read_afreq(afreq_path = afreq_path) %>%
    select(ID, ALT_FREQS)

# GeNA GWAS summary statistics
read_summstats <- function(celltype, analysis_name, resolution) {
    file_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats.txt")
    sumstats <- fread(file_path, select = 1:10) # read in results without the beta's
    sumstats[, `:=`(
        P = fifelse(as.numeric(P) == 0, .Machine$double.xmin * .Machine$double.eps, as.numeric(P)),
        celltype = celltype,
        permuted = "Real genotypes"
    )]
    return(sumstats)
}

# combine summary statistics for all cell types
sumstats_all_ct <- celltypes %>%
    map(\(celltype) read_summstats(celltype = celltype, analysis_name = analysis_name, resolution = resolution)) %>%
    rbindlist(fill = TRUE)

# format for plotting
sumstats_all_ct[, celltype := factor(str_replace(str_replace(celltype, pattern = "Other", replacement = "HSPC"), pattern = "_", replacement = " "), levels = unique(tenk_color_pal$major_cell_type)), ]
# remove SNPs with below 0.05 MAF
sumstats_all_ct <- merge(sumstats_all_ct, afreq, by = "ID", all.x = TRUE) %>%
    .[ALT_FREQS >= 0.05 & ALT_FREQS <= 0.95, ]
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

#-----------------------------------------------
# QQ plot with permuted data faceted by cell type
#-----------------------------------------------

combined_qqplot <- sumstats_all_ct %>%
    ggplot(aes(x = expected, y = observed, colour = celltype)) +
    labs(
        x = log10Pe,
        y = log10Po
    ) +
    # geom_point() + # can try scattermost for faster?
    geom_scattermore(pointsize = 5.3, pixels = c(1024, 1024)) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype = 2, colour = "red") +
    theme_classic() +
    scale_colour_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    theme(aspect.ratio = 1)

combined_qqplot %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/combined_qq_{analysis_name}_MAF_0.05.png"),
        width = 6, height = 5
    )

#-----------------------------------------------
# QQ plot with permuted data faceted by cell type
#-----------------------------------------------

# Permuted GeNA results
read_summstats_perm <- function(celltype, analysis_name, resolution) {
    file_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}_perm/GeNA_sumstats.txt")
    sumstats <- fread(file_path, select = 1:10) # read in results without the beta's
    sumstats[, `:=`(
        P = fifelse(as.numeric(P) == 0, .Machine$double.xmin * .Machine$double.eps, as.numeric(P)),
        celltype = celltype,
        permuted = "Permuted genotypes"
    )]
    return(sumstats)
}

sumstats_all_perm <- celltypes %>%
    map(\(celltype) read_summstats_perm(celltype = celltype, analysis_name = analysis_name, resolution = resolution)) %>%
    rbindlist(fill = TRUE)

# format for plotting
sumstats_all_perm[, celltype := factor(str_replace(str_replace(celltype, pattern = "Other", replacement = "HSPC"), pattern = "_", replacement = " "), levels = unique(tenk_color_pal$major_cell_type)), ]

# remove SNPs with below 0.05 MAF
sumstats_all_perm <- merge(sumstats_all_perm, afreq, by = "ID", all.x = TRUE) %>%
    .[ALT_FREQS >= 0.05 & ALT_FREQS <= 0.95, ]
setorder(sumstats_all_perm, P)
sumstats_all_perm[
    ,
    observed := -log10(P),
]
# calculate -log10 expected p value for each cell type
sumstats_all_perm[
    ,
    expected := -log10(ppoints(.N)),
    by = celltype
]

# create combined df with permuted and non-permuted data
sumstats_all <- rbindlist(list(sumstats_all_ct, sumstats_all_perm))

sumstats_all <- sumstats_all %>% .[
    ,
    permuted := factor(permuted, levels = c("Real genotypes", "Permuted genotypes")),
]

# faceted qqplot
combined_qqplot_facet <- sumstats_all %>%
    ggplot(aes(x = expected, y = observed, colour = permuted)) +
    facet_wrap(~celltype, ncol = 3) +
    labs(
        x = log10Pe,
        y = log10Po
    ) +
    geom_scattermore(pointsize = 11.3, pixels = c(1024, 1024)) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype = 2, colour = "red") +
    theme_classic() +
    # scale_colour_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    theme(
        aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = rel(0.8)),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8)),
        axis.text = element_text(size = rel(0.8)),
        axis.title = element_text(size = rel(1.2))
    )

combined_qqplot_facet %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/combined_qqplot_facet_{analysis_name}_MAF_0.05.pdf"),
        width = 8, height = 6
    )


#-----------------------------------------------
# Number of cells in cell type vs number of significant snps
#-----------------------------------------------

# NOTE: Making a newer version of this just showing lead (LD-pruned) SNPs.

# https://www.nature.com/articles/s41588-024-01909-1#:~:text=Therefore%2C%20we%20consider%20loci%20with%20GeNA%20P%E2%80%89%3C%E2%80%895%E2%80%89%C3%97%E2%80%8910%E2%88%928%20associations%20genome%2Dwide%20significant.
# GeNA manuscript: "we consider loci with GeNA P < 5 × 10−8 associations genome-wide significant."

sumstats_all_ct_sig_summary <- sumstats_all_ct[P < 5e-8, .N, by = celltype]

metadata <- get_latest_metadata()

ncells_nsnps_per_ct <- metadata %>%
    group_by(major_cell_type) %>%
    count() %>%
    rename("celltype" = major_cell_type, "Number of cells in cell type" = `n`) %>%
    left_join(sumstats_all_ct_sig_summary, by = "celltype") %>%
    rename("Number of significant csaQTLs" = N) %>%
    mutate(celltype = factor(celltype, levels = unique(tenk_color_pal$major_cell_type)))

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
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/combined_ncells_by_nsnps_plot_{analysis_name}.png"),
        width = 6, height = 5
    )

ncells_nsnps_per_ct %>%
    write_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/{resolution}/ncells_nsnps_per_ct_{analysis_name}.csv"))

# subset to just one SNP per LD-block? i.e. just the lead SNPs?

# -------------------------------
# save list of significant SNPs
# -------------------------------

sumstats_all_ct_sig <- sumstats_all_ct[P < 5e-8, ] %>%
    select(-observed, -expected)
sumstats_all_ct_sig %>%
    write_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/{resolution}/sumstats_all_ct_sig_{analysis_name}.csv"))

# -------------------------------
#  Copy of the above qq plot but with MAF 0.01 to demonstrate poor calibration at this MAF
# -------------------------------

# combine summary statistics for all cell types
MAF_0.01_sumstats_all_ct <- celltypes %>%
    map(\(celltype) read_summstats(celltype = celltype, analysis_name = analysis_name, resolution = resolution)) %>%
    rbindlist(fill = TRUE)

# format for plotting
MAF_0.01_sumstats_all_ct[, celltype := factor(str_replace(str_replace(celltype, pattern = "Other", replacement = "HSPC"), pattern = "_", replacement = " "), levels = unique(tenk_color_pal$major_cell_type)), ]
# remove SNPs with below 0.05 MAF
MAF_0.01_sumstats_all_ct <- merge(MAF_0.01_sumstats_all_ct, afreq, by = "ID", all.x = TRUE)
setorder(MAF_0.01_sumstats_all_ct, P)
MAF_0.01_sumstats_all_ct[
    ,
    observed := -log10(P),
]
# calculate -log10 expected p value for each cell type
MAF_0.01_sumstats_all_ct[
    ,
    expected := -log10(ppoints(.N)),
    by = celltype
]

# Permuted GeNA results
# combine summary statistics for all cell types
MAF_0.01_sumstats_all_ct <- celltypes %>%
    map(\(celltype) read_summstats(celltype = celltype, analysis_name = analysis_name, resolution = resolution)) %>%
    rbindlist(fill = TRUE)

# format for plotting
MAF_0.01_sumstats_all_ct[, celltype := factor(str_replace(str_replace(celltype, pattern = "Other", replacement = "HSPC"), pattern = "_", replacement = " "), levels = unique(tenk_color_pal$major_cell_type)), ]
# remove SNPs with below 0.05 MAF
MAF_0.01_sumstats_all_ct <- merge(MAF_0.01_sumstats_all_ct, afreq, by = "ID", all.x = TRUE)

setorder(MAF_0.01_sumstats_all_ct, P)

MAF_0.01_sumstats_all_ct[
    ,
    observed := -log10(P),
]
# calculate -log10 expected p value for each cell type
MAF_0.01_sumstats_all_ct[
    ,
    expected := -log10(ppoints(.N)),
    by = celltype
]
MAF_0.01_sumstats_all_perm <- celltypes %>%
    map(\(celltype) read_summstats_perm(celltype = celltype, analysis_name = analysis_name, resolution = resolution)) %>%
    rbindlist(fill = TRUE)

# format for plotting
MAF_0.01_sumstats_all_perm[, celltype := factor(str_replace(str_replace(celltype, pattern = "Other", replacement = "HSPC"), pattern = "_", replacement = " "), levels = unique(tenk_color_pal$major_cell_type)), ]

# remove SNPs with below 0.05 MAF
MAF_0.01_sumstats_all_perm <- merge(MAF_0.01_sumstats_all_perm, afreq, by = "ID", all.x = TRUE)
setorder(MAF_0.01_sumstats_all_perm, P)
MAF_0.01_sumstats_all_perm[
    ,
    observed := -log10(P),
]
# calculate -log10 expected p value for each cell type
MAF_0.01_sumstats_all_perm[
    ,
    expected := -log10(ppoints(.N)),
    by = celltype
]

# create combined df with permuted and non-permuted data
MAF_0.01_sumstats_all <- rbindlist(list(MAF_0.01_sumstats_all_ct, MAF_0.01_sumstats_all_perm))

MAF_0.01_sumstats_all <- MAF_0.01_sumstats_all %>% .[
    ,
    permuted := factor(permuted, levels = c("Real genotypes", "Permuted genotypes")),
]

log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

# faceted qqplot
combined_qqplot_facet <- MAF_0.01_sumstats_all %>%
    ggplot(aes(x = expected, y = observed, colour = permuted)) +
    facet_wrap(~celltype, ncol = 3) +
    labs(
        x = log10Pe,
        y = log10Po
    ) +
    geom_scattermore(pointsize = 11.3, pixels = c(1024, 1024)) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype = 2, colour = "red") +
    theme_classic() +
    # scale_colour_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    theme(
        aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(size = rel(1.2)),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(1.2)),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.4))
    )

combined_qqplot_facet %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/combined_qqplot_facet_{analysis_name}_MAF_0.01.pdf"),
        width = 8, height = 6
    )
