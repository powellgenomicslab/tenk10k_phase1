library(tidyverse)
library(viridis)
library(scales)
library(glue)
library(scattermore)
library(ggbeeswarm)
library(patchwork)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

celltype <- "Monocyte"
variant <- "12:9953308:T:TG"
name <- "CLEC12A"

# ⚙️ Functions ----

#' convert number/dose encoded genotypes to ref/alt bases in the SNP
geno_to_letters <- function(ref, alt, geno) {
    if (geno == 0) {
        out <- glue("{ref}/{ref}")
    } else if (geno == 1) {
        out <- glue("{ref}/{alt}")
    } else if (geno == 2) {
        out <- glue("{alt}/{alt}")
    } else if (geno == -9 || is.na(geno)) {
        out <- "Missing genotype"
    } else {
        stop(glue("Unexpected genotype value: {geno}"))
    }
    return(out)
}


# 📚 Read in the data ----

# UMAP data
## ADD THIS BACK IN WHEN I've SAVED UMAP csvs
# umap <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/{resolution}/{celltype}_umap") %>%
#     rename("barcode" = `...1`, "UMAP 1" = UMAP1, "UMAP 2" = UMAP2) %>%
#     mutate(
#         celltype = case_match(
#             celltype,
#             "NK1" ~ "NK_mature",
#             "NK3" ~ "NK_adaptive",
#             "NKint" ~ "NK_intermediate",
#             .default = celltype
#         )
#     )
csaQTL_npheno <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/major_cell_types/no_expr_pc_covars/{celltype}_neighbourhood_pheno.csv")) %>%
    rename("barcode" = `...1`)

# plot_data <- umap %>%
#     left_join(csaQTL_npheno, by = "barcode")

plot_data <- csaQTL_npheno

# boxplot data
csaQTL_spheno <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/major_cell_types/no_expr_pc_covars/{celltype}_sample_pheno.csv"))

# 🎻 violin box plots ----

# genotype boxplot for the locus - showing the sample-level phenotypes

# calculate the percentage of all major celltype cells belonging to each scpred subtype, per individual

ct_percentages <- plot_data %>%
    group_by(id) %>%
    mutate(total_cell_count_indiv = n()) %>%
    group_by(id, wg2_scpred_prediction) %>%
    mutate(
        ct_count_indiv = n(),
        pct_ct_indiv = ct_count_indiv / total_cell_count_indiv
    ) %>%
    select(id, wg2_scpred_prediction, pct_ct_indiv) %>%
    distinct() %>%
    pivot_wider(names_from = wg2_scpred_prediction, values_from = pct_ct_indiv, names_prefix = "pct_")

boxplot_data <- csaQTL_spheno %>%
    left_join(ct_percentages, by = "id")


ref <- str_split(variant, pattern = ":")[[1]][3]
alt <- str_split(variant, pattern = ":")[[1]][4]

sample_pheno <- boxplot_data %>%
    select(id, glue("spheno_{variant}"), {{ variant }}, sex, age, geno_PC1, geno_PC2, geno_PC3, geno_PC4, geno_PC5, geno_PC6, geno_PC7, BioHEART, starts_with("pct_")) %>%
    select(-pct_counts_mt) %>%
    mutate(across(starts_with("pct_"), ~ replace(., is.na(.), 0)))

# sanity check
# sample_pheno %>%
#     select(starts_with("pct_")) %>%
#     rowSums()

sample_pheno$Genotype <- sample_pheno %>%
    pull(variant) %>%
    map(\(geno) {
        geno_to_letters(ref, alt, geno)
    }) %>%
    unlist() %>%
    factor(levels = c(glue("{ref}/{ref}"), glue("{ref}/{alt}"), glue("{alt}/{alt}")))

ct_col <- tenk_color_pal %>%
    filter(major_cell_type == celltype) %>%
    select(color_major_cell_type) %>%
    distinct() %>%
    pull()

vlnplot_spheno <- ggplot(sample_pheno, aes(x = Genotype, y = !!sym(glue("spheno_{variant}")))) +
    geom_violin(color = NA, fill = ct_col, alpha = 0.7) +
    geom_boxplot(fill = NA, col = ct_col, width = 0.1, outlier.shape = NA, alpha = 1) +
    theme_classic() +
    theme(
        aspect.ratio = 1,
        text = element_text(size = 20),
        title = element_text(hjust = 0.5),
    ) +
    labs(y = "Sample-level phenotype", title = variant)

vlnplot_spheno %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_{name}_spheno_boxplot.pdf"),
        width = 5,
        height = 5
    )

# plot the raw percentages ----

pct_boxplot_data <- sample_pheno %>%
    pivot_longer(starts_with("pct_"), names_to = "Cell subtype", values_to = "Fraction of cells", names_prefix = "pct_")

vlnplot_cell_percentages <- ggplot(pct_boxplot_data, aes(x = Genotype, y = `Fraction of cells`)) +
    geom_violin(color = NA, fill = ct_col, alpha = 0.7) +
    geom_boxplot(fill = NA, col = ct_col, width = 0.1, outlier.shape = NA, alpha = 1) +
    facet_wrap(~`Cell subtype`, ncol = 2, scales = "free_y") +
    theme_classic() +
    theme(
        aspect.ratio = 1,
        text = element_text(size = 20),
        title = element_text(hjust = 0.5)
    )

vlnplot_cell_percentages %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_{name}_cell_percentages_boxplot.pdf"),
        width = 7.5,
        height = 7.5
    )

# plot residualised percentages, corrected for covariates ----

ct_pct_resid_list <- list()
for (ct in sample_pheno %>%
    select(starts_with("pct_")) %>%
    colnames()) {
    fit <- lm(as.formula(paste0(ct, "~ sex + age + geno_PC1 + geno_PC2 + geno_PC3 + geno_PC4 + geno_PC5 + geno_PC6 + geno_PC7 + BioHEART")), data = sample_pheno)
    ct_pct_resid_list[[paste0(ct, "_resid")]] <- residuals(fit) %>%
        as.data.frame()
}

ct_pct_resid_df <- ct_pct_resid_list %>%
    bind_cols()

names(ct_pct_resid_df) <- names(ct_pct_resid_list)

sample_pheno_resids <- bind_cols(sample_pheno, as_tibble(ct_pct_resid_df))

pct_boxplot_data_resids <- sample_pheno_resids %>%
    pivot_longer(starts_with("pct_"), names_to = "Cell subtype", values_to = "Fraction of cells", names_prefix = "pct_")

vlnplot_cell_percentages_resids <- ggplot(pct_boxplot_data_resids, aes(x = Genotype, y = `Fraction of cells`)) +
    geom_violin(color = NA, fill = ct_col, alpha = 0.7) +
    geom_boxplot(fill = NA, col = ct_col, width = 0.1, outlier.shape = NA, alpha = 1) +
    geom_quasirandom() +
    facet_wrap(~`Cell subtype`, ncol = 2, scales = "free_y") +
    theme_classic() +
    theme(
        aspect.ratio = 1,
        text = element_text(size = 20),
        title = element_text(hjust = 0.5)
    ) +
    labs(y = "Residualised fraction of cells")

vlnplot_cell_percentages_resids %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_{name}_cell_percentages_vln_residuals.pdf"),
        width = 14,
        height = 14
    )

# Individual box plots ---

boxplot_list <- list()

boxplot_list[["vlnplot_spheno"]] <- vlnplot_spheno

cell_subtypes <- pct_boxplot_data_resids$`Cell subtype` %>% unique()

# Add in the celltypes here to plot - use the ones with the greatest shifts in abundance
for (cell_type_plot in cell_subtypes[cell_subtypes %>% str_detect("_resid")]) {
    # filter data
    ct_pct_boxplot_data_resids <- pct_boxplot_data_resids %>%
        filter(`Cell subtype` == cell_type_plot)

    minor_ct <- cell_type_plot %>%
        str_remove("_resid") %>%
        str_replace_all("_", " ")

    minor_ct_col <- tenk_color_pal %>%
        filter(cell_type == minor_ct) %>%
        select(color) %>%
        distinct() %>%
        pull()

    # calculate boxplot stats manually to set axis limits
    axis_limits <- ct_pct_boxplot_data_resids %>%
        group_by(Genotype) %>%
        summarise(
            lower = quantile(`Fraction of cells`, 0.25) - 1.5 * IQR(`Fraction of cells`),
            upper = quantile(`Fraction of cells`, 0.75) + 1.5 * IQR(`Fraction of cells`)
        ) %>%
        ungroup() %>%
        summarise(
            lower = min(lower),
            upper = max(upper)
        )

    # coefs <- coef(lm(`Fraction of cells` ~ `15:39687137:C:T`, data = ct_pct_boxplot_data_resids))

    # don't plot the outliers as these skew the plot scales
    boxplot_cell_percentages_resids <- ggplot(ct_pct_boxplot_data_resids, aes(x = !!sym(variant), y = `Fraction of cells`, group = !!sym(variant))) +
        geom_violin(color = NA, fill = minor_ct_col, alpha = 0.7) +
        geom_boxplot(fill = NA, col = minor_ct_col, width = 0.1, outlier.shape = NA, alpha = 1) +
        # geom_beeswarm() +
        geom_smooth(method = "lm", col = minor_ct_col, aes(group = 1), se = FALSE) +
        theme_classic() +
        geom_smooth(method = "lm", se = TRUE) +
        theme(
            aspect.ratio = 1,
            text = element_text(size = 20),
            title = element_text(hjust = 0.5)
        ) +
        # scale_y_continuous(limits = c(axis_limits$lower, axis_limits$upper)) +
        scale_x_continuous(breaks = c(0, 1, 2), labels = levels(ct_pct_boxplot_data_resids$Genotype)) +
        labs(y = "Residualised fraction of cells", title = minor_ct)

    boxplot_list[[cell_type_plot]] <- boxplot_cell_percentages_resids
}

pct_boxplots <- patchwork::wrap_plots(boxplot_list, ncol = 3)

pct_boxplots %>%
    ggsave(
        filename = glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures/major_cell_types/umap/{celltype}_{name}_cell_percentages_vln_residuals_combined.pdf"),
        width = 12,
        height = 10
    )


