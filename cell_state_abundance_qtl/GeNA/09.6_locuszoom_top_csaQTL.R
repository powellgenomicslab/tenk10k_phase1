# micromamba activate r-python-dev

library(data.table)
library(tidyverse)
library(AnnotationHub)
library(locuszoomr)
library(glue)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")
an <- AnnotationHub(cache = "/home/blabow/.cache/R/AnnotationHub")
ensDb_v110 <- an[["AH113665"]]

# ----
analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

# celltype <- args[1]

celltype <- "NK"
# celltype <- "Monocyte"
# celltype <- "B"
# celltype <- "Dendritic"
# celltype <- "CD4_T"
# celltype <- "CD8_T" NO matching
# celltype <- "Unconventional_T"
# ----

# read in csaQTL raw p-values, these are the values we are plotting
# gena_summstats <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats.txt"))
afreq_path <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15.afreq"

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
    .[ALT_FREQS >= 0.05, ] # filter out low MAF results

plot_csaqtl_locuszoom <- function(variant, sumstats_csaqtl, celltype, tenk_color_pal) {
    plot_out_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/locuszoomR/top_csaQTL/{celltype}_{analysis_name}")

    locus_plot_list <- list() # make a list of plots to patchwork together, one csaQTL and one eqtl for each minor cell type

    # Get the SNP chr and pos
    # chr <- sumstats_csaqtl %>%
    #     filter(ID == {{ variant }}) %>%
    #     pull(CHROM)
    # pos <- sumstats_csaqtl %>%
    #     filter(ID == {{ variant }}) %>%
    #     pull(POS)
    print(variant)
    chr <- str_split_i(variant, pattern = ":", 1) %>%
        as.numeric()

    position <- str_split_i(variant, pattern = ":", 2) %>%
        as.numeric()

    ct_col <- tenk_color_pal %>%
        mutate(major_cell_type = str_replace(major_cell_type, " ", "_")) %>% # note remove later to format names nicely
        dplyr::filter(major_cell_type == {{ celltype }}) %>%
        dplyr::select(major_cell_type, color_major_cell_type) %>%
        dplyr::distinct() %>%
        dplyr::pull(color_major_cell_type)

    loc_csaqtl <- locus(
        data = sumstats_csaqtl,
        seqname = chr,
        xrange = c(position - 1e5, position + 1e5),
        # flank = 1e5,
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

    # add in the gene track
    locus_plot_list[[glue("gene_track")]] <- gg_genetracks(loc_csaqtl,
        gene_col = "gray",
        exon_col = "gray",
        exon_border = "gray",
    )

    combined_plots <- wrap_plots(locus_plot_list, ncol = 1)
    # combined_plots %>%
    #     ggsave(filename = glue("{plot_out_path}_locuszoom_caqtl_{variant}.png"), dpi = 600, width = 10, height = length(combined_plots) * 4)
    if (nchar(variant) > 30) {
        variant <- variant %>% str_trunc(30) # Truncate the really long variant names
    }
    combined_plots %>%
        ggsave(filename = glue("{plot_out_path}_locuszoom_caqtl_{variant}.pdf"), width = 5, height = length(combined_plots) * 2)

    return(combined_plots)
}


# get the lead SNPs
sumstats_gena_lead_snps <- fread(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats_lead_snps_MAF_0.05.tsv"))

sumstats_gena_lead_snps %>%
    slice_min(order_by = P, n = 10) %>%
    pull(ID) %>%
    map(\(variant) {
        plot_csaqtl_locuszoom(variant, sumstats_gena_ct, celltype, tenk_color_pal)
    })

# potential case examples ----

# Monocytes:
# plot_csaqtl_locuszoom("15:39687137:C:T", sumstats_gena_ct, celltype, tenk_color_pal)

# Dendritic cells:
# plot_csaqtl_locuszoom("7:119907178:ATG:A", sumstats_gena_ct, celltype, tenk_color_pal)
