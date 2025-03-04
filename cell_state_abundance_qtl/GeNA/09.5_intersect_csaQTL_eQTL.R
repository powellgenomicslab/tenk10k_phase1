library(data.table)
library(tidyverse)
library(glue)
library(rtracklayer)
library(GenomicRanges)
library(plyranges)
library(glue)
library(qvalue)

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

# read in csaQTL sig SNPs

# gena_sumstats_sig <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/combined/{resolution}/sumstats_all_ct_sig_{analysis_name}.csv"))
# gena_sumstats_sig <- gena_sumstats_sig %>%
#     mutate(CHR_POS = paste0(`#CHROM`, ":", POS))
celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt")

read_summstats <- function(celltype, analysis_name, resolution) {
    file_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats_lead_snps_MAF_0.05.tsv")
    sumstats <- fread(file_path, select = 1:10) # read in results without the beta's
    sumstats[, `:=`(
        P = fifelse(as.numeric(P) == 0, .Machine$double.xmin * .Machine$double.eps, as.numeric(P)),
        celltype = celltype
    )]
    return(sumstats)
}

# read in csaQTL lead snps
sumstats_lead_snps <- celltypes %>%
    map(\(celltype) read_summstats(celltype = celltype, analysis_name = analysis_name, resolution = resolution)) %>%
    rbindlist(fill = TRUE)

sumstats_gr <- makeGRangesFromDataFrame(
    df = sumstats_lead_snps,
    start.field = "POS",
    end.field = "POS",
    seqnames.field = "#CHROM",
    keep.extra.columns = TRUE
)
seqlevelsStyle(sumstats_gr) <- "UCSC"

# list of minor cell types to map functions across
minor_cell_types <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/unique_cell_types_wg2_scpred.txt")

# read in lead eQTLs, apply multiple testing correction per cell type
read_lead_eqtl <- function(minor_celltype) {
    eqtl_path <- glue("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/{minor_celltype}/{minor_celltype}_all_cis_cv_gene_level_results.tsv")

    eqtl_dt <- fread(eqtl_path,
        col.names = c("gene", "ACAT_p", "top_MarkerID", "top_pval"),
        colClasses = list(
            character = c("gene", "top_MarkerID"),
            numeric = c("ACAT_p", "top_pval")
        )
    )[
        , minor_celltype := minor_celltype # add column for cell type
    ][
        , qvalue := qvalue(ACAT_p)$qvalues # calculate q values
    ][
        qvalue < 0.05 # filter out non significant results
    ]

    return(eqtl_dt)
}

lead_eqtl <- minor_cell_types %>%
    map(\(minor_celltype) read_lead_eqtl(minor_celltype = minor_celltype)) %>%
    rbindlist()

# calculate largest nominal pval that is significant for each cell type
ct_pval_thresholds <- lead_eqtl[, max(top_pval), by = minor_celltype]
names(ct_pval_thresholds)[2] <- "max_pval_ct"

ct_pval_thresholds %>% write_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/eqtl_raw_pval_sig_thresholds.csv")

# read in all eqtls, filter all SNPs using those pvals
read_eqtl <- function(this_minor_celltype, ct_pval_thresholds) {
    ct_pval_thresh <- ct_pval_thresholds[minor_celltype == this_minor_celltype, max_pval_ct]
    eqtl_path <- glue("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/{this_minor_celltype}/{this_minor_celltype}_common_all_cis_raw_pvalues.tsv")
    eqtl_dt <- fread(
        eqtl_path,
        select = 2:19,
        col.names = c(
            "CHR", "POS", "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2",
            "MissingRate", "BETA", "SE", "Tstat", "var", "p.value", "p.value.NA",
            "Is.SPA", "N", "gene", "is_snp"
        ),
        colClasses = list(
            character = c("CHR", "MarkerID", "Allele1", "Allele2", "gene"),
            integer = c("POS", "N"),
            numeric = c("AC_Allele2", "AF_Allele2", "MissingRate", "BETA", "SE", "Tstat", "var", "p.value", "p.value.NA"),
            logical = "Is.SPA"
        )
    )[
        p.value < ct_pval_thresh # filter non-significant results using threshold based on the significant q-values above
    ][
        ,
        p.value := fifelse(p.value == 0, .Machine$double.xmin * .Machine$double.eps, p.value)
    ][
        ,
        "minor_celltype" := this_minor_celltype
    ]
    return(eqtl_dt)
}

all_eqtl_fdr_sig <- minor_cell_types %>%
    map(\(minor_celltype) {
        read_eqtl(this_minor_celltype = minor_celltype, ct_pval_thresholds = ct_pval_thresholds)
    }) %>%
    rbindlist()

# all_eqtl_fdr_sig[, max(p.value), by = minor_celltype]

# add in gene names

gencode_path <- "/directflow/SCCGGroupShare/projects/anncuo/reference_data/gencode.v44.basic.annotation.gtf"
gencode <- rtracklayer::import(gencode_path)
gencode_genes <- gencode[gencode$type == "gene", c("gene_id", "gene_name")] # filter down to just genes
gencode_genes$gene_id_noversion <- gencode_genes$gene_id %>% str_remove("\\.[0-9]+") # remove version numbers to match eqtl results
# intersect(gencode_genes$gene_id_noversion, all_eqtl_fdr_sig_gr$gene)
gene_name_map <- gencode_genes %>%
    as.data.frame() %>%
    select(gene_id_noversion, gene_name)

all_eqtl_fdr_sig <- all_eqtl_fdr_sig %>%
    left_join(gene_name_map, by = c("gene" = "gene_id_noversion"))

# convert to granges for intersection
all_eqtl_fdr_sig_gr <- makeGRangesFromDataFrame(
    df = all_eqtl_fdr_sig,
    start.field = "POS",
    end.field = "POS",
    seqnames.field = "CHR",
    keep.extra.columns = TRUE
)
seqlevelsStyle(all_eqtl_fdr_sig_gr) <- "UCSC"

csaQTL_overlapping_eQTL <- join_overlap_inner(sumstats_gr, all_eqtl_fdr_sig_gr, suffix = c("_GeNA", "_SAIGE"))

csaQTL_overlapping_eQTL %>%
    as.data.frame() %>%
    arrange(P) %>%
    fwrite(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/{resolution}_{analysis_name}_csaQTL_intersecting_eQTL.csv"))


# TODO: calculate the proportion of total csaQTL that are intersecting with eQTL?

# candidate case examples
# filtered to get more cell type specific eQTLs
csaQTL_overlapping_eQTL %>%
    as.data.frame() %>%
    as_tibble() %>%
    select(ID, celltype, minor_celltype, gene_name) %>%
    add_count(gene_name, name = "eqtl_celltype_count") %>%
    filter(eqtl_celltype_count <= 4) %>%
    arrange(eqtl_celltype_count) %>%
    fwrite(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/{resolution}_{analysis_name}_csaQTL_intersecting_eQTL_summary.csv"))

# TODO:

# Scatter Plot csaQTL vs eQTL -log10 p-values within a given locus should be corellated
