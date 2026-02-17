library(tidyverse)
library(data.table)
library(glue)

celltype <- "Monocyte"
analysis_name <- "single_variant_permutations"
resolution <- "major_cell_types"
variant <- "15:39687137:C:T"

# get GeNA GWAS summary statistics for plotting
read_summstats_perm <- function(perm_n, celltype, analysis_name, resolution, variant) {
    file_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/single_variant_{variant}_perm{perm_n}/GeNA_sumstats.txt")
    sumstats <- fread(file_path, select = 9) # return only the p-value # read in results without the beta's
    sumstats[, `:=`(
        P = fifelse(as.numeric(P) == 0, .Machine$double.xmin * .Machine$double.eps, as.numeric(P))
    )]
    return(sumstats)
}

# combine GeNA summary statistics for all cell types
perm_pvals <- 1:1000 %>%
    map({
        \(perm_n) read_summstats_perm(perm_n, celltype = celltype, analysis_name = analysis_name, resolution = resolution, variant = variant)
    }) %>%
    rbindlist()

# calculate an empirical p-value
p_val_real <- 6.62367982684486e-09

perm_pvals %>%
    arrange(P)

# get the number of permuted p-values below the real p-value
n_significant_perm <- perm_pvals %>%
    filter(P < p_val_real) %>%
    nrow()

empirical_pval <- n_significant_perm / 1000
empirical_pval

min(perm_pvals)

plot <- perm_pvals %>%
    ggplot(aes(x = P), breaks = seq(0, 1, 0.05)) +
    geom_histogram(bins = 20) +
    theme_classic() +
    scale_y_continuous(expand = expansion(c(0, 0.05)))

plot %>%
    ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/permutations/pval_hist.png")
