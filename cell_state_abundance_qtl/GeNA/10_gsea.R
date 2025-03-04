library(fgsea)
library(msigdbr)
library(tidyverse)
library(data.table)
library(glue)

celltype <- "NK"
resolution <- "major_cell_types"
analysis_name <- "no_expr_pc_covars"

# Retrieve human H (hallmark) gene set
msigdbr_df <- msigdbr(species = "human", category = "H")

# fixing format to work with fgsea
pathwaysH <- split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)

# read in the gene correlations for csaQTLs

var_gene_correlations <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_vargene_corrs.csv"))


# run fgsea enrichment
fgseaRes <- fgsea(pathways = pathwaysH, ranks, ..)
