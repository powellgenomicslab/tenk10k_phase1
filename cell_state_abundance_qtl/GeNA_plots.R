# Plot results for GeNA GWAS summary statistics

library(tidyverse)
library(qqman)
library(data.table)

sumstats <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/GeNA_sumstats.txt")

# Pvalue histogram
pval_hist <- sumstats %>%
    as_tibble() %>%
    ggplot(aes(x = P)) +
    geom_histogram(bins = 100)

pval_hist %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/NK_pval_hist.png", width = 5, height = 5)

# min(sumstats$P[sumstats$P != 0])

min(sumstats[P != 0, P])

# dim(sumstats[P == 0, ]) # replace zero's with r lowest float value

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

# Significant results
# sumstats_sig <- sumstats[P < 5e-8, , ][order(P)] %>%
#     as_tibble()

sumstats_no_zero_P <- sumstats[P != 1e-40, , ]
dim(sumstats_no_zero_P)

png(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/NK_manhattan_no_zero_P.png", width = 2800, height = 2000, res = 600)
manhattan(sumstats_no_zero_P,
    chr = "#CHROM", bp = "POS", snp = "ID", p = "P",
    main = "NK cells csaQTL", cex = 0.05, cex.axis = 0.4
)
dev.off()

# QQ plot
png(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/NK_qq_no_zero_P.png", width = 2000, height = 2000, res = 600)
qq(sumstats_no_zero_P$P,
    main = "NK cells csaQTL", cex = 0.3, cex.axis = 0.6
)
dev.off()
