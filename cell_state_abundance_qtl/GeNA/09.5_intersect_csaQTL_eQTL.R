library(data.table)
library(tidyverse)
library(glue)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

# read in csaQTL lead SNPs

gena_sumstats_sig <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/{resolution}/sumstats_all_ct_sig_{analysis_name}.csv"))
gena_sumstats_sig <- gena_sumstats_sig %>%
    mutate(CHR_POS = paste0(`#CHROM`, ":", POS))

sumstats_gr <- makeGRangesFromDataFrame(
    df = gena_sumstats_sig,
    start.field = "POS",
    end.field = "POS",
    seqnames.field = "#CHROM",
    keep.extra.columns = TRUE
)
seqlevelsStyle(sumstats_gr) <- "UCSC"

# read in all eQTL results
