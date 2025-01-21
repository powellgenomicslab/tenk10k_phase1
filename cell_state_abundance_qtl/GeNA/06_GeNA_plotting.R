# micromamba activate r-python-dev
library(tidyverse)
library(qqman)
library(data.table)
library(glue)

# ----
# Plot results for GeNA GWAS summary statistics
# ----

args <- commandArgs(trailingOnly = TRUE)
args <- list("CD4_TEM", "with_expr_pc_covars")

celltype <- args[1]
analysis_name <- args[2]

sumstats <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{celltype}/{analysis_name}/GeNA_sumstats.txt")) %>%
    .[, P := as.numeric(P), ] %>%
    .[, P := ifelse(P == 0, .Machine$double.xmin * .Machine$double.eps, P)]

# swap values reported as 0 for the smallest possible floating point.
# NOTE: numbers below .Machine$double.xmin (2.225074e-308) are reported with reduced precision
# for max precision I should report anything below 2.225074e-308 as P < 2.225074e-308?
# For now I will keep values below this, noting that precision for these values is reduced

# min(sumstats$P[sumstats$P != 0])
# min(sumstats$P)

# Pvalue histogram
pval_hist <- sumstats %>%
    ggplot(aes(x = P)) +
    geom_histogram(breaks = seq(0, 1, 0.05)) +
    scale_y_continuous(expand = expansion(c(0, 0.05)))

pval_hist %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{celltype}_{analysis_name}_pval_hist.png"),
        width = 5, height = 5
    )

# min(sumstats[,P])

# manhattan plot
png(
    filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{celltype}_{analysis_name}_manhattan.png"),
    width = 2800, height = 2000, res = 600
)
manhattan(sumstats,
    chr = "#CHROM", bp = "POS", snp = "ID", p = "P",
    main = glue("{celltype} csaQTL"), cex = 0.05, cex.axis = 0.4
)
dev.off()

# QQ plot
png(
    filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{celltype}_{analysis_name}_qq.png"),
    width = 2000, height = 2000, res = 600
)
qq(sumstats$P,
    main = glue("{celltype} csaQTL"), cex = 0.3, cex.axis = 0.6
)
dev.off()

# results passing genome-wide significance threshold
sumstats_sig <- sumstats[P < 5e-8, , ][order(P)]

# save them
sumstats_sig %>%
    fwrite(
        glue(
            "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{celltype}/{analysis_name}/GeNA_sumstats_sig.csv"
        )
    )


# QTL box plots
