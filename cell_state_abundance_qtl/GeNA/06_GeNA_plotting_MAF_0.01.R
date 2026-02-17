# micromamba activate r-python-dev
library(tidyverse)
library(data.table)
library(glue)
library(scattermore)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

# ----
# Plot results for GeNA GWAS summary statistics
# ----

args <- commandArgs(trailingOnly = TRUE)
# args <- list("major_cell_types", "NK", "no_expr_pc_covars")

res <- args[1]
celltype <- args[2]
analysis_name <- args[3]

plot_out_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{res}/{celltype}_{analysis_name}_")

# swap values reported as 0 for the smallest possible floating point.
# NOTE: numbers below .Machine$double.xmin (2.225074e-308) are reported with reduced precision
# for max precision I should report anything below 2.225074e-308 as P < 2.225074e-308?
# For now I will keep values below this, noting that precision for these values is reduced

sumstats <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{res}/{celltype}/{analysis_name}/GeNA_sumstats.txt")) %>%
    .[, P := as.numeric(P), ] %>%
    .[, P := ifelse(P == 0, .Machine$double.xmin * .Machine$double.eps, P)] %>%
    .[, permuted := "Real genotypes", ]

sumstats_permuted <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{res}/{celltype}/{analysis_name}_perm/GeNA_sumstats.txt")) %>%
    .[, P := as.numeric(P), ] %>%
    .[, P := ifelse(P == 0, .Machine$double.xmin * .Machine$double.eps, P), ] %>%
    .[, permuted := "Permuted genotypes", ]

sumstats_combined <- sumstats %>%
    bind_rows(sumstats_permuted)

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

# NOTE: permutation results indicated that calibration was poor for previous threshold of 1%

sumstats <- sumstats %>%
    left_join(afreq, by = "ID")

sumstats_permuted <- sumstats_permuted %>%
    left_join(afreq, by = "ID")

sumstats_combined <- sumstats_combined %>%
    left_join(afreq, by = "ID")

# min(sumstats$P[sumstats$P != 0])
# sumstats %>% filter(P == min(sumstats$P))
# sumstats_permuted %>% filter(P == min(sumstats_permuted$P))


# Pvalue histogram
pval_hist <- sumstats %>%
    ggplot() +
    geom_histogram(data = sumstats, aes(x = P, fill = permuted), breaks = seq(0, 1, 0.05), alpha = 0.4) +
    geom_histogram(data = sumstats_permuted, aes(x = P, fill = permuted), breaks = seq(0, 1, 0.05), alpha = 0.4) +
    scale_y_continuous(expand = expansion(c(0, 0.05))) +
    theme_classic() +
    labs(title = glue("{celltype} {analysis_name}")) +
    theme(aspect.ratio = 1)

pval_hist %>%
    ggsave(
        filename = glue("{plot_out_path}_pval_hist_MAF_0.01.png"),
        width = 5, height = 5
    )

setorder(sumstats_combined, P)
sumstats_combined[
    ,
    observed := -log10(P),
]
# calculate -log10 expected p value for each cell type
sumstats_combined[
    ,
    expected := -log10(ppoints(.N)),
    by = permuted
]

log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

qqplot <- sumstats_combined %>%
    ggplot(aes(x = expected, y = observed, colour = permuted), size = 3.4) +
    labs(
        x = log10Pe,
        y = log10Po
    ) +
    labs(title = glue("{celltype} {analysis_name}")) +
    # geom_point() + # can try scattermost for faster?
    geom_scattermore(pointsize = 5.3, pixels = c(1024, 1024)) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype = 2, colour = "red") +
    theme_classic() +
    theme(aspect.ratio = 1)
# scale_colour_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$wg2_scpred_prediction))

qqplot %>%
    ggsave(
        filename = glue("{plot_out_path}_qqplot_MAF_0.01.png"),
        width = 5, height = 5
    )


# check correllation between permuted and real
# table(sumstats$ID == sumstats_permuted$ID)
# cor(sumstats$P, sumstats_permuted$P, method = "spearman")
# permuted and real are not correllated?

# calculate p-value lambda inflation factor
# chisq <- qchisq(1 - sumstats$P, 1)
# lambda <- median(chisq) / qchisq(0.5, 1)
# lambda
# sumstats_combined_sig <- sumstats_combined[P < 5e-8, , ][order(P)]
# table(sumstats_combined_sig$permuted)
