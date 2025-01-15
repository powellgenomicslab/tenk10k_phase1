library(data.table)
library(dplyr)
library(ggplot2)

# All eQTL results from the SAIGE-QTL run (December 2024 freeze)
# The data is organised by cell type
saige_dir = '/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/'
celltypes = list.files(saige_dir)
celltypes = celltypes[!(celltypes %in% c('README.md','CD4_TCM_sample_perm0','annotate','coloc_results'))]

# cell type colours from Blake
cols_df = as.data.frame(fread('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/colour_palette_table.tsv'))

# Manhattan plot function
manhattan <- function(df0, col="#94B9D4"){
    # enforce numeric position
    df0$POS = as.numeric(df0$POS)

    # prepare dataframe for Manhattan
    don <- df0 %>% 

      # Compute chromosome size
      group_by(CHR) %>% 
      summarise(chr_len=max(POS)) %>% 

      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      select(-chr_len) %>%

      # Add this info to the initial dataset
      left_join(df0, ., by=c("CHR"="CHR")) %>%

      # Add a cumulative position of each SNP
      arrange(CHR, POS) %>%
      mutate(BPcum=POS+tot)

    # axes
    axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
    p = ggplot(don, aes(x=BPcum, y=-log10(p.value))) +
    
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("#D0CECE", col), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks=axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    xlab('Genomic position') + 
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      text = element_text(size=20),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    )
    return(p)
}

options(repr.plot.width = 12, repr.plot.height = 6) 
for (celltype in celltypes){
    col = cols_df[cols_df$wg2_scpred_prediction == celltype,'color']
#     print(col)
    pvals_file = paste0(saige_dir,celltype,'/',celltype,'_common_all_cis_raw_pvalues.tsv')
    pvals_df = fread(pvals_file)
    p = manhattan(pvals_df, col=col)
    print(p+ggtitle(celltype))
}
