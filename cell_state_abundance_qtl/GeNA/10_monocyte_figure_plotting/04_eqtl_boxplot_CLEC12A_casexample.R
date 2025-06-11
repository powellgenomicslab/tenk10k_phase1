library(tidyverse)
library(viridis)
library(scales)
library(glue)
library(scattermore)
library(ggbeeswarm)
library(patchwork)
library(data.table)
source("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

# parameters
minor_celltype <- "CD16_Mono"
celltype <- "Monocyte"
gene <- "ENSG00000172322"
gene_name <- "CLEC12A"
snp <- "12:9953308:T:TG"

# boxplot data (genotype data)
csaQTL_spheno <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/major_cell_types/no_expr_pc_covars/{celltype}_sample_pheno.csv"))

# read in npheno for the barcode-individual mapping
metadata <- get_latest_metadata() %>%
    select(barcode, cpg_id, individual)

# read in data for the particular gene
boxplots_data_dir <- "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/celltype_specificity/"
genotypes_dir <- paste0(boxplots_data_dir, "boxplots_genotypes/")
expression_dir <- paste0(boxplots_data_dir, "boxplots_expression/")
expr_filename <- paste0(expression_dir, minor_celltype, "_", gene, ".csv")
# geno_filename <- paste0(genotypes_dir, gsub(":", "_", snp), ".csv")

df_e <- as.data.frame(fread(expr_filename))
colnames(df_e)[1] <- "cell"
# df_g <- as.data.frame(fread(geno_filename))
# df_g$V1 <- c()
df_e_inds <- inner_join(df_e, adata_obs[, c("cell", "individual")])
colnames(df_e_inds)[2] <- "gene"
# calculate mean expression
df_e_mean <- df_e_inds %>%
    group_by(individual) %>%
    summarize(gene_mean = mean(gene, na.rm = TRUE))
# combine genotype and mean expression data
df_to_plot <- inner_join(df_e_mean, df_g)
df_to_plot <- df_to_plot[rowSums(is.na(df_to_plot)) == 0, ]
colnames(df_to_plot)[3] <- "snp"


# create box plot
p <- ggplot(df_to_plot, aes(x = as.factor(snp), y = gene_mean))
p <- p + geom_violin(color = ct1_col, fill = ct1_col, alpha = 0.7)
p <- p + geom_boxplot(width = 0.1, outlier.shape = NA, col = ct1_col, fill = ct1_col)
p <- p + theme_classic() + theme(text = element_text(size = 20))
gene_name <- gene_info[gene_info$gene == gene, ]$gene_name
p <- p + ylab(paste0(gene_name, " mean")) + ggtitle(minor_celltype) + labs(title = gene_name, x = snp)

p %>% ggsave(filename = , width = 4, height = 5)
