# mm activate r-python-dev

library(data.table)
library(tidyverse)
library(glue)

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

# read in the curated genelist
deficiency_genes <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/intersections/csaQTL_gene_list_celltype_deficiency.csv") %>%
    distinct()
# covert gene symbols that don't match up with gencode names
# from online gene name converter tool
name_conversion <- tribble(
    ~Symbol, ~gene_name,
    "RECQL3", "BLM",
    "CD3Z", "CD247",
    "CD21", "CR2",
    "SLP76", "LCP2",
    "CD20", "MS4A1",
    "NBS1", "NBN",
    "POLE1", "POLE",
    "SKIV2L", "SKIC2",
    "TTC37", "SKIC3",
    "ERBB21P", "ERBIN"
)
deficiency_genes <- deficiency_genes %>%
    left_join(name_conversion, by = "Symbol") %>%
    mutate(gene_name = fifelse(is.na(gene_name), Symbol, gene_name))

# read in the genic csaQTL snps
snps_intersecting_gene_bodies <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/genic_snps/{resolution}_{analysis_name}_snps_intersecting_gene_bodies.csv"))
snps_within_100kb_cis_window <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/genic_snps/{resolution}_{analysis_name}_snps_within_100kb_cis_window_of_genes.csv"))

# find csaQTL intersecting the deficiency genes
deficiency_gene_intersecting_csaqtl <- snps_intersecting_gene_bodies %>%
    inner_join(
        deficiency_genes,
        by = "gene_name"
    ) %>%
    mutate(SNP_within_gene_body = TRUE) %>%
    mutate(SNP_within_100kb_cis_window = TRUE)

deficiency_gene_intersecting_csaqtl %>%
    select(gene_name, Phenotype, celltype, ID, P) %>%
    distinct()

# find csaQTL within 100kb of deficiency genes
deficiency_genes_in_100kb_csaqtl <- snps_within_100kb_cis_window %>%
    inner_join(
        deficiency_genes,
        by = "gene_name"
    ) %>%
    mutate(SNP_within_100kb_cis_window = TRUE)

deficiency_genes_in_100kb_csaqtl %>%
select(gene_name, Phenotype, celltype, ID, P) %>%
    distinct()

all_intersecting_summary <- deficiency_gene_intersecting_csaqtl %>%
    bind_rows(deficiency_genes_in_100kb_csaqtl) %>%
    distinct(gene_name, ID, celltype, .keep_all = TRUE) %>%
    arrange(gene_name, P) %>%
    mutate(SNP_within_gene_body = if_else(is.na(SNP_within_gene_body), FALSE, SNP_within_gene_body))

all_intersecting_summary %>%
    write_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/{resolution}_{analysis_name}_csaQTL_intersecting_curated_deficiency_genes.csv"))

# Note: some of the provided gene names may be different to the gencode annotation
## Check how many of the provided symbols are not in the annotation

# library(rtracklayer)
# library(GenomicRanges)

# get gene where names are different
# gencode_path <- "/directflow/SCCGGroupShare/projects/anncuo/reference_data/gencode.v44.basic.annotation.gtf"
# gencode <- rtracklayer::import(gencode_path)
# gencode_genes <- gencode[gencode$type == "gene", ]
# setdiff(unique(deficiency_genes$Symbol), gencode$gene_name)
# setdiff(unique(deficiency_genes$gene_name), gencode$gene_name)

# setdiff(unique(name_conversion$gene_name), gencode_genes$gene_name)

#  [1] "CD3Z"
#  [2] "SLP76"
#  [3] "RECQL3"
#  [4] "NBS1"
#  [5] "POLE1"
#  [6] "11q23del"
#  [7] "Del10p13-p14"
#  [8] "Large (3Mb) deletion of 22q11.2"
#  [9] "Unknown / environment"
# [10] "ERBB21P"
# [11] "SKIV2L"
# [12] "TTC37"
# [13] " Unknown"
# [14] "CD20"
# [15] "CD21"
# [16] "  Unknown"
# [17] "Mutation or chromosomal deletion at 14q32"
