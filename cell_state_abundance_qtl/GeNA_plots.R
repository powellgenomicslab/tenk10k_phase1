# micromamba activate r-python-dev
# Plot results for GeNA GWAS summary statistics

library(tidyverse)
library(qqman)
library(data.table)

sumstats <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/GeNA_sumstats.txt")

# Pvalue histogram
pval_hist <- sumstats %>%
    as_tibble() %>%
    ggplot(aes(x = P)) +
    geom_histogram(breaks = seq(0, 1, 0.05)) +
    scale_y_continuous(expand = expansion(c(0, 0.05)))

pval_hist %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/NK_pval_hist.png", width = 5, height = 5)

# min(sumstats[,P])

# manhattan plot
png(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/NK_manhattan.png", width = 2800, height = 2000, res = 600)
manhattan(sumstats,
    chr = "#CHROM", bp = "POS", snp = "ID", p = "P",
    main = "NK cells csaQTL", cex = 0.05, cex.axis = 0.4
)
dev.off()

# QQ plot
png(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/NK_qq.png", width = 2000, height = 2000, res = 600)
qq(sumstats$P,
    main = "NK cells csaQTL", cex = 0.3, cex.axis = 0.6
)
dev.off()

# results passing genome-wide significance threshold

sumstats_sig <- sumstats[P < 5e-8, , ][order(P)] %>%
    as_tibble()
