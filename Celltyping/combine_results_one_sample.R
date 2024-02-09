library(data.table)
library(ggplot2)
library(Seurat)
library(viridis)

args <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

scpred_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/"
samples = list.files(scpred_dir, pattern = "S.")
sample = samples[i]

out_file = paste0(scpred_dir, sample, "/combined_metadata.csv")
# if (file.exists(out_file) == TRUE){quit(save = "no", status = 1, runLast = FALSE)}

print(paste0("Combine cell classification results for ",sample))

azimuth_file = paste0(scpred_dir, sample, "/step2_azimuth/azimuth.RDS")
if (file.exists(azimuth_file) == FALSE){quit()}
azimuth_obj = readRDS(azimuth_file)
# load hierarchical scPred output file
hierscpred_file = paste0(scpred_dir, sample, "/step3_hierscpred/hier_scpred.RDS")
if (file.exists(hierscpred_file) == FALSE){quit()}
hierscpred_obj = readRDS(hierscpred_file)
# extract metadata
df1 = azimuth_obj@meta.data
df2 = hierscpred_obj@meta.data
# combine and save
df3 = cbind(df1,df2)
write.csv(df3, out_file)

x = df3[,c("predicted.celltype.l2", "scpred_prediction")]

createheatmap <- function(x, order = TRUE){

  xaxis <- x[, 1]
  yaxis <- x[, 2]

  cont_table <- table(xaxis, yaxis)
  prop_table <- cont_table / rowSums(cont_table)
  cont_table <-  data.frame(cont_table)
  prop_table <-  data.frame(prop_table)
  cont_table$Prop <- prop_table$Freq

  axes_names <- names(x)
  grey_col <-  grey.colors(5000, start = 1, end = 0) # Color of the text in the heatmap
  grey_col[1:2500] <- grey_col[1]
  grey_col[2501:5000] <- grey_col[5000]


  if(order){
    label_order <- c("CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "Treg", "CD4 Proliferating",
                     "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating", "gdT", "MAIT", "ILC",
                     "dnT", "NK", "NK_CD56bright", "NK Proliferating", "B naive", "B intermediate",
                     "B memory", "Plasmablast", "CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC",
                     "ASDC", "HSPC", "Platelet", "Eryth", "Doublet")

    cont_table[,1] <- factor(cont_table[,1], levels = label_order)
    cont_table[,2] <- factor(cont_table[,2], levels = label_order)
  }
  # Colormap
  hmcol <- viridis(500)
  grey_col[1] <- hmcol[1] #Plot 0 as invisible

  p <- ggplot(cont_table, aes(x = xaxis, y = yaxis, fill = Prop)) +
    geom_tile(colour = "white", size = 0.25) +
    geom_text(aes(label = Freq, colour = Prop), size = 4) +
    scale_colour_gradientn(colours = grey_col, limits = c(0, 1.1), guide = "none") +
    scale_fill_gradientn(colours = hmcol,
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         labels = c(0, 0.25, 0.5, 0.75, 1),
                         guide = "none") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.background = element_blank(),
          panel.border = element_blank(),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    xlab(axes_names[1]) +
    ylab(axes_names[2])

  p2 <- ggplot(cont_table, aes(x = xaxis, y = yaxis, fill = Prop)) +
    geom_tile(colour = "white", size = 0.25) +
    geom_text(aes(label = round(Prop, 2), colour = Prop), size = 4) +
    scale_colour_gradientn(colours = grey_col, limits = c(0, 1.1), guide = "none") +
    scale_fill_gradientn(colours = hmcol,
                         limits = c(0, 1),
                         breaks = c(0, 0.25, 0.5, 0.75, 1),
                         labels = c(0, 0.25, 0.5, 0.75, 1),
                         guide = "none") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.background = element_blank(),
          panel.border = element_blank(),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    xlab(axes_names[1]) +
    ylab(axes_names[2])

  list(count = p, prop = p2, cont_table = cont_table)

}    

fig_dir <- paste0(scpred_dir, "figures/")
pdf(paste0(fig_dir, sample"_count_comparison_heatmap.pdf"), width=20, height=20)
createheatmap(x)[1]
dev.off()

pdf(paste0(fig_dir, sample"_prop_comparison_heatmap.pdf"), width=20, height=20)
createheatmap(x)[2]
dev.off()

