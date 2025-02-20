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


# read in the Gene positions
gencode_path <- "/directflow/SCCGGroupShare/projects/anncuo/reference_data/gencode.v44.basic.annotation.gtf"
gencode <- rtracklayer::import(gencode_path)

# filter down to just gene bodies (introns + exons)
gencode_genes <- gencode[gencode$type == "gene", ]
# get csaQTL SNPs that intersect with the gene bodies

snps_intersecting_gene_bodies <- join_overlap_left(sumstats_gr, gencode_genes)

snps_intersecting_gene_bodies %>%
    as.data.frame() %>%
    mutate(intersecting_gene_body = if_else(is.na(gene_name), TRUE, FALSE)) %>% 
    write_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/genic_snps/{resolution}_{analysis_name}_snps_intersecting_gene_bodies.csv"))

snps_within_100kb_cis_window <- join_overlap_left(sumstats_gr, gencode_genes, maxgap = 100000)

snps_within_100kb_cis_window %>%
    as.data.frame() %>%
    arrange(P)
    mutate(intersecting_gene_body_100kb_window = if_else(!is.na(gene_name), TRUE, FALSE)) %>% 
    write_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/genic_snps/{resolution}_{analysis_name}_snps_within_100kb_cis_window_of_genes.csv"))

# 96% of csaQTL SNPs are within 100kb of gene bodies
# length(unique(snps_within_100kb_cis_window$ID)) / length(unique(sumstats_gr$ID))
snps_intersecting_gene_bodies %>%
    as.data.frame() %>%
    mutate(intersecting_gene_body = if_else(!is.na(gene_name), TRUE, FALSE)) %>%
    filter(intersecting_gene_body) %>%
    group_by(celltype) %>%
    mutate(
        n_snps_intersecting_genes = n_distinct(ID)
    ) %>%
    select(celltype, n_snps_intersecting_genes) %>%
    distinct() %>% 
    arrange(desc(n_snps_intersecting_genes))

snps_intersecting_gene_bodies %>%
    as.data.frame() %>%
    mutate(intersecting_gene_body = if_else(!is.na(gene_name), TRUE, FALSE)) %>%
    filter(intersecting_gene_body == TRUE) %>%
    arrange(P) %>%
    select(celltype, gene_name, P, ALT_FREQS, ID) %>%
    group_by(gene_name) %>%
    slice_min(P) %>% 
    arrange(P) %>% print(n=20)
