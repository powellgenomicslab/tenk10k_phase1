# micromamba activate r-python-dev

library(data.table)
library(tidyverse)
library(AnnotationHub)
library(locuszoomr)
library(glue)
library(patchwork)

source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")
an <- AnnotationHub(cache = "/home/blabow/.cache/R/AnnotationHub")
ensDb_v110 <- an[["AH113665"]]

# ----
analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"
# celltype <- "NK"
celltype <- "Monocyte"
# celltype <- "B"
# celltype <- "Dendritic"
# celltype <- "CD4_T"
# celltype <- "CD8_T" NO matching
# celltype <- "Unconventional_T"
# ----

# read in the top csaQTL that are intersecting with eQTL, these are the genes / loci we are interested in plotting
sumstats_sig_csaQTL_eQTL <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/{resolution}_{analysis_name}_csaQTL_intersecting_eQTL.csv"))

sumstats_sig_csaQTL_eQTL_ct <- sumstats_sig_csaQTL_eQTL %>%
    dplyr::filter(celltype == {{ celltype }})

# read in csaQTL raw p-values, these are the values we are plotting
gena_summstats <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats.txt"))

read_afreq <- function(afreq_path) {
    afreq_df <- readr::read_tsv(
        afreq_path,
        comment = "#",
        col_names = c("CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT"),
        col_types = cols(
            CHROM = col_character(),
            ID = col_character(),
            REF = col_character(),
            ALT = col_character(),
            ALT_FREQS = col_double(),
            OBS_CT = col_integer(),
        )
    )
    return(afreq_df)
}

afreq_path <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15.afreq"
afreq <- read_afreq(afreq_path = afreq_path) %>%
    dplyr::select(ID, ALT_FREQS)

# get GeNA GWAS summary statistics for plotting
read_summstats <- function(celltype, analysis_name, resolution) {
    file_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats.txt")
    sumstats <- fread(file_path, select = 1:10) # read in results without the beta's
    sumstats[, `:=`(
        P = fifelse(as.numeric(P) == 0, .Machine$double.xmin * .Machine$double.eps, as.numeric(P)),
        celltype = celltype
    )]
    return(sumstats)
}

# combine GeNA summary statistics for all cell types
sumstats_gena_ct <- celltype %>%
    read_summstats(analysis_name = analysis_name, resolution = resolution)

sumstats_gena_ct <- merge(sumstats_gena_ct, afreq, by = "ID", all.x = TRUE) %>%
    .[ALT_FREQS >= 0.05 & ALT_FREQS <= 0.95, ] # filter out low MAF results

# get minor celltypes corresponding to the specified major cell type we are plotting
major_minor_ct_map <- tenk_color_pal %>%
    mutate(major_cell_type = str_replace(major_cell_type, " ", "_")) %>%
    dplyr::select(major_cell_type, wg2_scpred_prediction) %>%
    dplyr::filter(major_cell_type == celltype)

# get eqtl summary statistics, we are plotting these p values
read_eqtl_raw_pvals <- function(this_minor_celltype) {
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
        ,
        p.value := fifelse(p.value == 0, .Machine$double.xmin * .Machine$double.eps, p.value) # handle small values
    ][
        ,
        "minor_celltype" := this_minor_celltype # column specifying minor cell type
    ]
    return(eqtl_dt)
}
# read in saige summary statistics for all of the minor celltypes corresponding to the major cell type we are looking at
sumstats_saige_ct <- major_minor_ct_map$wg2_scpred_prediction %>%
    map(\(minor_celltype) {
        read_eqtl_raw_pvals(this_minor_celltype = minor_celltype)
    }) %>%
    rbindlist()

# add gene names to saige output
gencode_path <- "/directflow/SCCGGroupShare/projects/anncuo/reference_data/gencode.v44.basic.annotation.gtf"
gencode <- rtracklayer::import(gencode_path)
gencode_genes <- gencode[gencode$type == "gene", c("gene_id", "gene_name")] # filter down to just genes
gencode_genes$gene_id_noversion <- gencode_genes$gene_id %>% str_remove("\\.[0-9]+") # remove version numbers to match eqtl results
# intersect(gencode_genes$gene_id_noversion, all_eqtl_fdr_sig_gr$gene)
gene_name_map <- gencode_genes %>%
    as.data.frame() %>%
    dplyr::select(gene_id_noversion, gene_name)

sumstats_saige_ct <- sumstats_saige_ct %>%
    left_join(gene_name_map, by = c("gene" = "gene_id_noversion"))

ct_pval_thresholds <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/eqtl_raw_pval_sig_thresholds.csv")

plot_gene_locuszoom <- function(gene, sumstats_csaqtl, sumstats_eqtl, celltype, minor_cell_types = NULL, ct_pval_thresholds, tenk_color_pal) {
    plot_out_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/locuszoomR/{celltype}_{analysis_name}")

    locus_plot_list <- list() # make a list of plots to patchwork together, one csaQTL and one eqtl for each minor cell type

    # log10P <- expression(paste("-log"[10], plain(P))) # for axis labels

    ct_col <- tenk_color_pal %>%
        mutate(major_cell_type = str_replace(major_cell_type, " ", "_")) %>% # note remove later to format names nicely
        dplyr::filter(major_cell_type == {{ celltype }}) %>%
        dplyr::select(major_cell_type, color_major_cell_type) %>%
        dplyr::distinct() %>%
        dplyr::pull(color_major_cell_type)

    loc_csaqtl <- locus(
        data = sumstats_csaqtl,
        gene = gene,
        flank = 1e5,
        ens_db = ensDb_v110,
        p = "P",
        labs = "ID",
        pos = "POS"
    )
    locus_plot_list[[glue("{celltype} csaQTL")]] <- gg_scatter(
        loc_csaqtl,
        labels = "index",
        ylab = paste0(" GeNA ", "-log10(P) "),
        scheme = c("grey", ct_col, "red"),
        xticks = FALSE
    ) + annotate("text", x = -Inf, y = 5e-8, col = "grey", label = glue("  P = 5e-8\n", .trim = FALSE), vjust = 0, hjust = 0) +
        annotate("text", x = Inf, y = Inf, col = "black", label = glue("{str_replace(celltype, '_', ' ')}  "), vjust = 1, hjust = 1)

    for (minor_celltype in minor_cell_types) {
        eqtl_plot_data <- sumstats_eqtl %>%
            dplyr::filter(minor_celltype == {{ minor_celltype }}) %>%
            dplyr::filter(gene_name == {{ gene }})

        # get significance threshold used in the cell type
        ct_pval_thresh <- ct_pval_thresholds %>%
            dplyr::filter(minor_celltype == {{ minor_celltype }}) %>%
            dplyr::pull(max_pval_ct)

        ct_col_minor <- tenk_color_pal %>%
            dplyr::filter(wg2_scpred_prediction == minor_celltype) %>% # note swap cols later for nicer formatting
            dplyr::pull(color)

        if (nrow(eqtl_plot_data) > 0) { # don't plot if no eqtl results for the cell type / gene combination (not all genes are in the eqtl results)
            loc_eqtl <- locus(
                data = eqtl_plot_data,
                gene = gene,
                flank = 1e5,
                ens_db = ensDb_v110,
                p = "p.value",
                labs = "MarkerID",
                pos = "POS"
            )

            locus_plot_list[[minor_celltype]] <- gg_scatter(loc_eqtl,
                pcutoff = ct_pval_thresh,
                labels = "index",
                scheme = c("grey", ct_col_minor, "red"),
                ylab = paste0(gene, " SAIGE ", "-log10(P) "),
                xticks = FALSE
            ) +
                annotate("text", x = -Inf, y = -log10(ct_pval_thresh), col = "grey", label = glue("  P = {round(ct_pval_thresh, 4)}\n", .trim = FALSE), vjust = 0, hjust = 0) +
                annotate("text", x = Inf, y = Inf, col = "black", label = glue("{str_replace(minor_celltype, '_', ' ')}  "), vjust = 1, hjust = 1)
        }
    }
    locus_plot_list[[glue("gene_track")]] <- gg_genetracks(loc_csaqtl,
        highlight = gene,
        highlight_col = "cornflowerblue",
        gene_col = "gray",
        exon_col = "gray",
        exon_border = "gray",
    )

    combined_plots <- wrap_plots(locus_plot_list, ncol = 1)
    combined_plots %>%
        ggsave(filename = glue("{plot_out_path}_{gene}_locuszoom_caqtl_eqtl.png"), dpi = 600, width = 10, height = length(combined_plots) * 4)
    ggsave(filename = glue("{plot_out_path}_{gene}_locuszoom_caqtl_eqtl.pdf"), width = 5, height = length(combined_plots) * 2)
}

# TODO: also plot eQTL in a different panel
intersecting_matching_celltypes <- sumstats_sig_csaQTL_eQTL_ct %>%
    # just get genes that are eqtl in the same cell type
    dplyr::filter(minor_celltype %in% major_minor_ct_map$wg2_scpred_prediction) %>%
    dplyr::arrange(P)

intersecting_matching_celltypes %>%
    write_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/intersections/{resolution}_{analysis_name}_{celltype}_csaQTL_intersecting_eQTL_matching_celltypes.csv"))

# get top ten significant csaQTL that are eqtl genes
genes_plot <- intersecting_matching_celltypes %>%
    dplyr::group_by(gene_name) %>%
    dplyr::summarize(P = min(P)) %>%
    dplyr::arrange(P) %>%
    # dplyr::slice_min(P, n = 10) %>%
    dplyr::pull(gene_name)

# sumstats_sig_csaQTL_eQTL_ct %>%
#     group_by(gene_name, minor_celltype) %>%
#     summarize(P = min(P)) %>%
#     arrange(P)
# NOTE: currenlty doesn't work for genes

genes_plot[!str_detect(genes_plot, "ENSG")] %>% map(\(gene) {
    plot_gene_locuszoom(gene, sumstats_gena_ct, sumstats_saige_ct, celltype, major_minor_ct_map$wg2_scpred_prediction, ct_pval_thresholds, tenk_color_pal)
})



# genes_plot_monocytes <- c(
#     "CDC42BPA",
#     "SPRED2",
#     "MBP",
#     "RXRA",
#     "AOAH",
#     "LYZ"
# )
# # LYZ antisense is also significant csaQTL and eQTL
# # can't plot it currently because it's not in the annotation

# genes_plot_monocytes %>% map(\(gene) {
#     plot_gene_locuszoom(gene, sumstats_gena_ct, sumstats_saige_ct, celltype, major_minor_ct_map$wg2_scpred_prediction, ct_pval_thresholds)
# })

# genes_plot_bcells <- c(
#     "IGHM",
#     "IGHG3"
# )

# genes_plot_bcells %>% map(\(gene) {
#     plot_gene_locuszoom(gene, sumstats_gena_ct, sumstats_saige_ct, celltype, major_minor_ct_map$wg2_scpred_prediction)
# })


# plot_gene_locuszoom("KLF2", sumstats_gena_ct, sumstats_saige_ct, celltype, minor_cell_types = major_minor_ct_map$wg2_scpred_prediction)


#################################################################

# plot_gene_locuszoom <- function(gene, sumstats_ct, celltype) {
#     loc_csaqtl <- locus(data = sumstats_ct, gene = gene, flank = 1e5, ens_db = ensDb_v110, p = "P", labs = "ID", pos = "POS")

#     plot_out_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/locuszoomR/{celltype}_{analysis_name}")

#     png(
#         filename = glue("{plot_out_path}_{gene}_locuszoom.png"),
#         width = 3200, height = 2800, res = 600
#     )

#     locus_plot(loc_csaqtl,
#         highlight = gene,
#         gene_col = "gray",
#         exon_col = "gray",
#         exon_border = "gray",
#         highlight_col = "cornflowerblue",
#         ylab = paste0("-log10(P) ", celltype, " GeNA ", gene),
#     )

#     dev.off()
# }
