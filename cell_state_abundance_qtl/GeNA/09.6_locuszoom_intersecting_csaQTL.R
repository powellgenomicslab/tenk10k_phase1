library(data.table)
library(tidyverse)
library(AnnotationHub)
library(locuszoomr)

an <- AnnotationHub(cache = "/home/blabow/.cache/R/AnnotationHub")
ensDb_v110 <- an[["AH113665"]]

gene <- "KLF2"
celltype <- "NK"

gena_summstats <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/major_cell_types/NK/no_expr_pc_covars/GeNA_sumstats.csv")
gena_summstats <- gena_summstats[ALT_FREQS <= 0.05, ]


loc <- locus(data = df_to_plot, gene = gene_name, flank = 1e5, ens_db = ensDb_v110)
