# micromamba activate r-python-dev

library(data.table)
library(tidyverse)

# I think the ID must be old pos from GrCH37??
gena_tenk10k_NK <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/GeNA_sumstats.txt") %>%
    mutate(chr_pos_hg38 = paste0(`#CHROM`, ":", `POS`))


# 1. insect my results with the onek1k GeNA results ----

# significant loci published in rumker et al from the onek1k
# seems to be grch37-based
gena_onek1k <- tribble(
    ~`Cell Type`, ~CHR, ~POS, ~rsID, ~`Effect Allele`, ~`Other Allele`, ~`Cyto. Band`, ~MAF, ~P,
    "Myeloid", 15, 80263217, "rs3826007", "T", "C", "15q25.1", 0.27, 2.61e-08,
    "NK", 2, 111851212, "rs13025330", "T", "C", "2q13", 0.22, 1.76e-09,
    "NK", 11, 128070535, "rs519062", "G", "A", "11q24.3", 0.28, 2.48e-13,
    "NK", 12, 10583611, "rs3003", "T", "C", "12p13.2", 0.18, 1.96e-11,
    "NK", 19, 16441973, "rs56133626", "A", "G", "19p13.11", 0.33, 1.96e-13
) %>% mutate(ID = paste0(`CHR`, ":", `POS`))

intersecting <- gena_tenk10k_NK %>% inner_join(gena_onek1k, by = "ID", suffix = c("_tenk10k", "_onek1k"))
intersecting %>% select(`#CHROM`, `POS_tenk10k`, P_tenk10k, `Cell Type`)
intersecting %>% colnames()

# 2. Identify eQTL that are also csaQTL (simple overlap) ----

# read in eQTL results
eqtl_tenk1k_NK <- fread("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/associatr/NK_gene_tests_raw_pvals.txt")
colnames(eqtl_tenk1k_NK) <- c("CHR", "POS", "gene", "raw_p_val")

# get lead SNPs
eqtl_lean_snp_only <- eqtl_tenk1k_NK[, .SD[which.min(raw_p_val)], by = gene]

# get FDR - corrected data
eqtl_fdr <- fread("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/associatr/acat_fdr_gene_files/NK.tsv")

# join to FDR-coorected eQTL, filter for significant eqtl only
eqtl_sig <- eqtl_lean_snp_only %>%
    left_join(eqtl_fdr, by = "gene") %>%
    .[qv < 0.05, , ] %>%
    mutate(chr_pos_hg38 = paste0(str_remove(CHR, "chr"), ":", `POS`))

# get the significant GeNA associations
gena_tenk10k_sig <- gena_tenk10k_NK[P < 5e-8, , ]

candidate_csa_e_qtl <- gena_tenk10k_sig %>%
    inner_join(eqtl_sig, by = "chr_pos_hg38")

candidate_csa_e_qtl %>% select(chr_pos_hg38, gene, P, qv)
