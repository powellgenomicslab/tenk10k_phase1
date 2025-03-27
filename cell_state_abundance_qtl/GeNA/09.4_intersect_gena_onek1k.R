library(data.table)
library(tidyverse)
library(glue)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)
library(janitor)

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

# Intersection with all csaQTL SNPs (Lead + LD SNPs) ----

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
    as.data.frame() %>%
    select(rsID, ID_tenk10k, P_tenk10k, celltype, `Cell.Type`, ID_onek1k)


overlapping_gap <- join_overlap_inner(sumstats_gr, gena_sumstats_sig_onek1k, maxgap = 1000, suffix = c("_tenk10k", "_onek1k"))

# overlapping_gap %>%
#     as.data.frame()

# get lead snps corresponding to each onek1k locus - use these for plotting n_phenos
overlapping_gap %>%
    as.data.frame() %>%
    # janitor::clean_names() %>%
    group_by(celltype, `Cell.Type`, rsID) %>%
    slice_min(P_tenk10k) %>%
    select(rsID, ID_tenk10k, P_tenk10k, celltype, `Cell.Type`, ID_onek1k)

# Intersection with lead snps -----

summstats_lead_snps <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/major_cell_types/summstats_all_lead_csaqtl.csv")

summstats_lead_snps_gr <- makeGRangesFromDataFrame(
    df = summstats_lead_snps,
    start.field = "POS",
    end.field = "POS",
    seqnames.field = "#CHROM",
    keep.extra.columns = TRUE
)
seqlevelsStyle(summstats_lead_snps_gr) <- "UCSC"

# # no direct overlap with lead SNPs
# overlapping_lead_snps <- join_overlap_inner(summstats_lead_snps_gr, gena_sumstats_sig_onek1k, suffix = c("_tenk10k", "_onek1k"))

# these are the lead SNPs that I plot in the supplementary figure:
overlapping_lead_snps <- join_overlap_inner(summstats_lead_snps_gr, gena_sumstats_sig_onek1k, maxgap = 100000, suffix = c("_tenk10k", "_onek1k"))

overlapping_lead_snps %>%
    as.data.frame() %>%
    # janitor::clean_names() %>%
    group_by(celltype, `Cell.Type`, rsID) %>%
    slice_min(P_tenk10k) %>%
    select(rsID, CHR_POS_hg38, ID_tenk10k, P_tenk10k, celltype, `Cell.Type`, ID_onek1k)


# overlapping_gap %>%
#     as_tibble() %>%
#     filter(rsID == "rs13025330") %>%
#     arrange(celltype, ID_tenk10k) %>%
#     select(celltype, `Cell.Type`, ID_tenk10k, P_tenk10k, CHR_POS_hg38) %>%
#     print(n = 100)


##############################################################################################################
# LD calculations for manuscript // replication supp tables

# # calculate LD between the onek1k csaQTL and corresponding lead variants from tenk10k analysis:
# plink2 --pfile /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15 --ld 11:128200640:A:G 11:128246346:T:TA
# # result : r^2=0.942694

# plink2 --pfile /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15 --ld 15:79970875:C:T 15:79972101:G:A
# # result: r^2=0.991216

# plink2 --pfile /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15 --ld 19:16331162:G:A 19:16331208:G:A
# # result: r^2=0.995559

# # calculate LD for the replicated (non-intersecting) onek1k csaQTL and closest lead variant from my analysis:
# plink2 --pfile /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15 --ld 2:111093635:C:T 2:111050100:G:A
# # result" 0.455808

###################################################################################################################
