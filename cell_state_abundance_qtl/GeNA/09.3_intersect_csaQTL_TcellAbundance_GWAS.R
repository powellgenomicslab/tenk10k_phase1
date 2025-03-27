# mm activate r-python-dev

library(data.table)
library(tidyverse)
library(glue)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)

# parameters ------------------------
analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"
# -----------------------------------

# ----------------------------------------------------------
# Identify SNPS within gene body or cis window around genes
# ----------------------------------------------------------

# 📚 read in the data ----

# read in summary statistics for the significant GeNA results

gena_sumstats_sig <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/{resolution}/sumstats_all_ct_sig_{analysis_name}.csv"))
# convert to granges

sumstats_gr <- makeGRangesFromDataFrame(
    df = gena_sumstats_sig,
    start.field = "POS",
    end.field = "POS",
    seqnames.field = "#CHROM",
    keep.extra.columns = TRUE
)
seqlevelsStyle(sumstats_gr) <- "UCSC"

# read in the T cell abundance GWAS summary statistics from https://doi.org/10.1038/s41467-024-51095-1,
# downloaded from https://doi.org/10.5281/zenodo.12582912

tca_meta_sumstats <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/intersections/tcell_fraction_gwas/TOPMed_AoU_Meta_analysis.txt")
tca_meta_sumstats_sig <- tca_meta_sumstats[P < 5e-8, ]

tca_meta_sumstats_sig_gr <- makeGRangesFromDataFrame(
    df = tca_meta_sumstats_sig,
    start.field = "BP",
    end.field = "BP",
    seqnames.field = "CHR",
    keep.extra.columns = TRUE
)

seqlevelsStyle(tca_meta_sumstats_sig_gr) <- "UCSC"

# Get the csaQTL overlapping the SNPs associated with T cell abundance phenotypes
csaQTL_overlapping_tca_GWAS_snps <- join_overlap_inner(sumstats_gr, tca_meta_sumstats_sig_gr, suffix = c("_GeNA", "_TOPMed_AoU_Meta_analysis"))

csaQTL_overlapping_tca_GWAS_snps %>%
    as_tibble() %>%
    select(celltype, ID, P_GeNA) %>%
    arrange(P_GeNA) %>%
    print(n = 100)

csaQTL_overlapping_tca_GWAS_snps %>%
    as.data.frame() %>%
    write_csv(
        glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/{resolution}_{analysis_name}_csaQTL_intersecting_tca_gwas.csv")
    )
