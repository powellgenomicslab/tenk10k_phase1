# micromamba activate r-python-dev
library(tidyverse)
library(qqman)
library(data.table)
library(glue)
library(scattermore)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

# ----
# Plot results for GeNA GWAS summary statistics
# ----

args <- commandArgs(trailingOnly = TRUE)
args <- list("major_cell_types", "Dendritic", "no_expr_pc_covars")

res <- args[1]
celltype <- args[2]
analysis_name <- args[3]

plot_out_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{res}/{celltype}_{analysis_name}_")

sumstats <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{res}/{celltype}/{analysis_name}/GeNA_sumstats.txt")) %>%
    .[, P := as.numeric(P), ] %>%
    .[, P := ifelse(P == 0, .Machine$double.xmin * .Machine$double.eps, P)] %>%
    .[, permuted := "Real genotypes", ]

sumstats_permuted <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{res}/{celltype}/{analysis_name}_perm/GeNA_sumstats.txt")) %>%
    .[, P := as.numeric(P), ] %>%
    .[, P := ifelse(P == 0, .Machine$double.xmin * .Machine$double.eps, P), ] %>%
    .[, permuted := "Permuted genotypes", ]

# check correllation between permuted and real
# table(sumstats$ID == sumstats_permuted$ID)
# cor(sumstats$P, sumstats_permuted$P, method = "spearman")
# permuted and real are not correllated

# calculate p-value lambda inflation factor
# chisq <- qchisq(1 - sumstats$P, 1)
# lambda <- median(chisq) / qchisq(0.5, 1)
# lambda

sumstats_combined <- sumstats %>%
    bind_rows(sumstats_permuted)

# swap values reported as 0 for the smallest possible floating point.
# NOTE: numbers below .Machine$double.xmin (2.225074e-308) are reported with reduced precision
# for max precision I should report anything below 2.225074e-308 as P < 2.225074e-308?
# For now I will keep values below this, noting that precision for these values is reduced

# min(sumstats$P[sumstats$P != 0])
sumstats %>% filter(P == min(sumstats$P))
sumstats_permuted %>% filter(P == min(sumstats_permuted$P))


# Pvalue histogram
pval_hist <- sumstats %>%
    ggplot() +
    geom_histogram(data = sumstats, aes(x = P, fill = permuted), breaks = seq(0, 1, 0.05), alpha = 0.4) +
    geom_histogram(data = sumstats_permuted, aes(x = P, fill = permuted), breaks = seq(0, 1, 0.05), alpha = 0.4) +
    scale_y_continuous(expand = expansion(c(0, 0.05))) +
    theme_classic() +
    labs(title = glue("{celltype} {analysis_name}"))

pval_hist %>%
    ggsave(
        filename = glue("{plot_out_path}_pval_hist.png"),
        width = 5, height = 5
    )

# min(sumstats[,P])

# manhattan plot
# png(
#     filename = glue("{plot_out_path}_manhattan.png"),
#     width = 2800, height = 2000, res = 600
# )
# manhattan(sumstats,
#     chr = "#CHROM", bp = "POS", snp = "ID", p = "P",
#     main = glue("{celltype} csaQTL"), cex = 0.05, cex.axis = 0.4
# )
# dev.off()

# QQ plot
# png(
#     filename = glue("{plot_out_path}_qq.png"),
#     width = 2000, height = 2000, res = 600
# )
# qq(sumstats$P,
#     main = glue("{celltype} csaQTL"), cex = 0.3, cex.axis = 0.6
# )
# dev.off()

# calculate observed and expected -log10P

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
    geom_scattermore() +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5, linetype = 2, colour = "red") +
    geom_hline(yintercept = -log10(5e-8)) +
    theme_classic() # +
# scale_colour_manual(values = setNames(tenk_color_pal$color, tenk_color_pal$wg2_scpred_prediction))

qqplot %>%
    ggsave(
        filename = glue("{plot_out_path}_qqplot.png"),
        width = 5, height = 5
    )

# results passing genome-wide significance threshold
sumstats_sig <- sumstats[P < 5e-8, , ][order(P)]
# sumstats_combined_sig <- sumstats_combined[P < 5e-8, , ][order(P)]
# table(sumstats_combined_sig$permuted)

# https://www.nature.com/articles/s41588-024-01909-1#:~:text=Therefore%2C%20we%20consider%20loci%20with%20GeNA%20P%E2%80%89%3C%E2%80%895%E2%80%89%C3%97%E2%80%8910%E2%88%928%20associations%20genome%2Dwide%20significant.
# GeNA manuscript: "we consider loci with GeNA P < 5 × 10−8 associations genome-wide significant."


# save them
sumstats_sig %>%
    fwrite(
        glue(
            "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{res}/{celltype}/{analysis_name}/GeNA_sumstats_sig.csv"
        )
    )

