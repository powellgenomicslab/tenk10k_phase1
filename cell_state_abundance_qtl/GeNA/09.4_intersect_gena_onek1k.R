library(data.table)
library(tidyverse)
library(glue)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

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

# data from https://doi.org/10.1038/s41588-024-01909-1
# Supplementary Table 1.
chain <- import.chain("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/intersections/hg19_liftover/hg19ToHg38.over.chain")
gena_sumstats_sig_onek1k <- tribble(
    ~`Cell Type`, ~CHR, ~POS, ~rsID, ~`Effect Allele`, ~`Other Allele`, ~`Cyto. Band`, ~MAF, ~P,
    "Myeloid", 15, 80263217, "rs3826007", "T", "C", "15q25.1", 0.27, 2.61e-08,
    "NK", 2, 111851212, "rs13025330", "T", "C", "2q13", 0.22, 1.76e-09,
    "NK", 11, 128070535, "rs519062", "G", "A", "11q24.3", 0.28, 2.48e-13,
    "NK", 12, 10583611, "rs3003", "T", "C", "12p13.2", 0.18, 1.96e-11,
    "NK", 19, 16441973, "rs56133626", "A", "G", "19p13.11", 0.33, 1.96e-13
) %>%
    mutate(ID = paste0(`CHR`, ":", `POS`)) %>%
    makeGRangesFromDataFrame(
        start.field = "POS",
        end.field = "POS",
        seqnames.field = "CHR",
        keep.extra.columns = TRUE
    )
seqlevelsStyle(gena_sumstats_sig_onek1k) <- "UCSC"

gena_sumstats_sig_onek1k <- liftOver(gena_sumstats_sig_onek1k, chain) %>%
    unlist()

gena_sumstats_sig_onek1k$CHR_POS_hg38 <- paste0(str_remove(seqnames(gena_sumstats_sig_onek1k)@values, "chr"), ":", ranges(gena_sumstats_sig_onek1k))

overlapping <- join_overlap_inner(sumstats_gr, gena_sumstats_sig_onek1k, suffix = c("_tenk10k", "_onek1k"))

overlapping %>%
    as.data.frame()

overlapping_gap <- join_overlap_inner(sumstats_gr, gena_sumstats_sig_onek1k, maxgap = 100000, suffix = c("_tenk10k", "_onek1k"))

# overlapping_gap %>%
#     as.data.frame()

overlapping_gap %>%
    as.data.frame() %>%
    janitor::clean_names() %>%
    group_by(celltype, `cell_type`, rs_id) %>%
    slice_min(p_tenk10k)
