library(tidyverse)
library(data.table)
library(glue)
library(scattermore)
library(ggnewscale)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

analysis_name <- "no_expr_pc_covars"
resolution <- "major_cell_types"

# 📚 read in the data ----

# celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/unique_cell_types_wg2_scpred.txt")
celltypes <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt")
celltypes <- celltypes[celltypes != "ALL"]

# Get the minor allele frequencies
read_afreq <- function(afreq_path) {
    read_tsv(
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
}
afreq_path <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15.afreq"
afreq <- read_afreq(afreq_path = afreq_path) %>%
    dplyr::select(ID, ALT_FREQS)

# GeNA GWAS summary statistics
read_summstats <- function(celltype, analysis_name, resolution) {
    file_path <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats.txt")
    sumstats <- fread(file_path, select = 1:10) # read in results without the beta's
    sumstats[, `:=`(
        P = fifelse(as.numeric(P) == 0, .Machine$double.xmin * .Machine$double.eps, as.numeric(P)),
        celltype = celltype
    )]
    return(sumstats)
}

# combine summary statistics for all cell types
sumstats_all_ct <- celltypes %>%
    map(\(celltype) read_summstats(celltype = celltype, analysis_name = analysis_name, resolution = resolution)) %>%
    rbindlist(fill = TRUE)

# format for plotting
sumstats_all_ct[, celltype := factor(str_replace(celltype, pattern = "_", replacement = " "), levels = unique(tenk_color_pal$major_cell_type)), ]
sumstats_all_ct <- merge(sumstats_all_ct, afreq, by = "ID", all.x = TRUE) %>%
    .[ALT_FREQS >= 0.05, ]

setorder(sumstats_all_ct, P)
sumstats_all_ct[
    ,
    neg_log10_P := -log10(P),
]

# downsample the non-significant points for faster plotting
# table(sumstats_all_ct$P < 0.0001)

sumstats_all_ct_sig <- sumstats_all_ct %>%
    dplyr::filter(P < 0.001)

# keep only 10% of data points with P below 0.001
sumstats_all_ct_notsig <- sumstats_all_ct %>%
    dplyr::filter(P >= 0.001) %>%
    group_by(`#CHROM`) %>%
    slice_sample(prop = 0.1)

sumstats_subset <- bind_rows(sumstats_all_ct_sig, sumstats_all_ct_notsig) %>%
    rename("CHROM" = `#CHROM`)

plot_data <- sumstats_subset %>%
    # Compute chromosome size
    group_by(CHROM) %>%
    summarise(chr_len = as.numeric(max(POS))) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(sumstats_subset, ., by = c("CHROM" = "CHROM")) %>%
    # Add a cumulative position of each SNP
    arrange(CHROM, POS) %>%
    mutate(BPcum = POS + tot)

axisdf <- plot_data %>%
    group_by(CHROM) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)

log10P <- expression(paste("-log"[10], plain(P)))

# test <- plot_data %>% slice_head(n=10000)

# for single combined plot
# width <- 8
# height <- 3
# px <- 350
# ptsize <- 9.2

# for split plot
# width <- 12
# height <- 14
px <- 1250
ptsize <- 7.7
aspect_ratio <- 0.25

all_ct_manhattan <- ggplot(data = plot_data, aes(x = BPcum, y = neg_log10_P)) +
    # Show all points
    geom_scattermore(
        data = plot_data %>% dplyr::filter(P >= 5e-8),
        aes(color = as.factor(CHROM)), pointsize = ptsize, pixels = c(px, px * aspect_ratio),
        show.legend = FALSE
    ) +
    # geom_point(aes(color = as.factor(CHROM)), size = 1) +
    scale_color_manual(values = rep(c("grey", "lightgrey"), 22)) +
    new_scale_color() +
    # highlight significant loci
    geom_scattermore(
        data = plot_data %>% dplyr::filter(P < 5e-8),
        aes(color = celltype), pointsize = ptsize, pixels = c(px, px * aspect_ratio)
    ) +
    scale_colour_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$center, expand = expansion(mult = .02), guide = guide_axis(check.overlap = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + # remove space between plot area and x axis
    geom_hline(yintercept = -log10(5e-8), linetype = 2) +
    labs(
        y = log10P,
        x = "Chromosome"
    ) +
    guides(color = guide_legend("Major cell type", override.aes = list(size = 6, shape = 15)), fill = "none") +
    # Custom the theme:
    theme_classic() +
    theme(
        # legend.position = "none",
        aspect.ratio = aspect_ratio,
        # panel.border = element_blank(),
        strip.background = element_rect(fill = NA, color = NA),
        # strip.background = element_blank(),
        # strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        element_blank()
    )


# all_ct_manhattan %>% ggsave(
#     filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/combined_manhattan.png"),
#     width = 8, height = 3
# )

# all_ct_manhattan %>% ggsave(
#     filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/combined_manhattan.pdf"),
#     width = 8, height = 3
# )


ct_split_manhattan <- all_ct_manhattan +
    facet_wrap(~celltype, ncol = 2)

ct_split_manhattan %>% ggsave(
    filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/split_combined_manhattan.pdf"),
    width = 10, height = 8,
)

# make new column that can be used to colour by chromosome and cell type
plot_data <- plot_data %>%
    mutate(ct_chr = paste0(celltype, "_", CHROM))

chr_pal <- tenk_color_pal %>%
    select(major_cell_type, color_major_cell_type_light, color_major_cell_type_light2) %>%
    # pivot_longer(cols = -major_cell_type, names_to = "dark_or_light", values_to = "chr_cols") %>%
    distinct()

# Make colour palette for non-significant loci, with alternating light-dark for chromosomes
colorvec <- c()
# gross but it works
for (ct in tenk_color_pal$major_cell_type) {
    for (chr in 1:22) {
        if (chr %% 2 == 0) {
            colorvec[glue("{ct}_{chr}")] <- chr_pal[chr_pal$major_cell_type == ct, "color_major_cell_type_light"]
        } else {
            colorvec[glue("{ct}_{chr}")] <- chr_pal[chr_pal$major_cell_type == ct, "color_major_cell_type_light2"]
        }
    }
}

all_ct_manhattan_alt_cols <- ggplot(data = plot_data, aes(x = BPcum, y = neg_log10_P)) +
    # Show all points
    geom_scattermore(
        data = plot_data %>% dplyr::filter(P >= 5e-8),
        aes(color = as.factor(ct_chr)), pointsize = ptsize, pixels = c(px, px * aspect_ratio),
        show.legend = FALSE
    ) +
    # geom_point(aes(color = as.factor(CHROM)), size = 1) +
    scale_color_manual(values = colorvec) +
    new_scale_color() +
    # highlight significant loci
    geom_scattermore(
        data = plot_data %>% dplyr::filter(P < 5e-8),
        aes(color = celltype), pointsize = ptsize, pixels = c(px, px * aspect_ratio)
    ) +
    scale_colour_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$center, expand = expansion(mult = .02), guide = guide_axis(check.overlap = TRUE)) +
    # scale_y_continuous(expand = expansion(mult = c(.1, .1))) + # remove space between plot area and x axis
    geom_hline(yintercept = -log10(5e-8), linetype = 2) +
    labs(
        y = log10P,
        x = "Chromosome"
    ) +
    guides(color = guide_legend("Major cell type", override.aes = list(size = 6, shape = 15)), fill = "none") +
    # Custom the theme:
    theme_classic() +
    theme(
        # legend.position = "none",
        aspect.ratio = aspect_ratio,
        # panel.border = element_blank(),
        strip.background = element_rect(fill = NA, color = NA),
        # strip.background = element_blank(),
        # strip.text = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        element_blank()
    ) +
    facet_wrap(~celltype, ncol = 2)

all_ct_manhattan_alt_cols %>% ggsave(
    filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/split_combined_manhattan_alternative_colours.pdf"),
    width = 10, height = 8,
)
