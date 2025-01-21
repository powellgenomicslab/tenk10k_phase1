library(tidyverse)
library(qqman)
library(data.table)
library(glue)


analysis_name <- "with_expr_pc_covars"


# 📚 read in the data ----

celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/unique_cell_types_wg2_scpred.txt")

read_summstats <- function(celltype, analysis_name) {
    sumstats <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{celltype}/{analysis_name}/GeNA_sumstats.txt")) %>%
        .[, P := as.numeric(P), ] %>%
        .[, P := ifelse(P == 0, .Machine$double.xmin * .Machine$double.eps, P)] %>%
        .[, celltype := celltype]
    return(sumstats)
}
# combine summary statistics for all cell types
sumstats_all_ct <- celltypes %>%
    map(\(celltype) read_summstats(celltype = celltype, analysis_name = analysis_name)) %>%
    rbindlist(fill = TRUE)

# ---------------------------
# qq plot for all cell types
# ---------------------------

sumstats_all_ct %>%
    ggplot2(aes())

#-----------------------------------------------
# Number of cells in cell type vs number of significant snps?
#-----------------------------------------------

# subset to just one SNP per LD-block? i.e. just the lead SNPs?
