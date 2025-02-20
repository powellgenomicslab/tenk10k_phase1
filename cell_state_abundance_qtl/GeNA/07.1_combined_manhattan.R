library(tidyverse)
library(qqman)
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
    select(ID, ALT_FREQS)

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
    filter(P < 0.001)

# keep only 10% of data points with P below 0.001
sumstats_all_ct_notsig <- sumstats_all_ct %>%
    filter(P >= 0.001) %>%
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
    select(-chr_len) %>%
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

width <- 10
height <- 6

all_ct_manhattan <- ggplot( data = plot_data, aes(x = BPcum, y = neg_log10_P)) +
    # Show all points
    geom_scattermore(
        data = plot_data %>% filter(P >= 5e-8),
        aes(color = as.factor(CHROM)), pointsize = 9.2, pixels = c(300 * width, 300 * height),
         show.legend = FALSE
    ) +
    # geom_point(aes(color = as.factor(CHROM)), size = 1) +
    scale_color_manual(values = rep(c("grey", "lightgrey"), 22)) +
    new_scale_color() +
    # highlight significant loci 
    geom_scattermore(
        data = plot_data %>% filter(P < 5e-8),
        aes(color = celltype), pointsize = 9.2, pixels = c(300 * width, 300 * height)
    ) +
    scale_colour_manual(values = setNames(unique(tenk_color_pal$color_major_cell_type), unique(tenk_color_pal$major_cell_type))) +
    # custom X axis:
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$center, expand = expansion(mult = .02)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + # remove space between plot area and x axis
    geom_hline(yintercept = -log10(5e-8), linetype = 2) +
    # Custom the theme:
    theme_classic() +
    theme(
        # legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    labs(
        y = log10P,
        x = "Chromosome"
    ) + 
    guides(color = guide_legend("Major cell type"), fill = "none")


all_ct_manhattan %>% ggsave(
    filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/{resolution}/combined_plots/combined_manhattan.png"),
     width = width, height = height
)

