library(data.table)
library(ggplot2)

cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
samples = list.files(cellranger_dir, pattern = "S.")

df_to_plot = data.frame()
for (sample in samples){
    cellranger_metrics_file1 = paste0(cellranger_dir, sample, "/outs/metrics_summary.csv")
    cellranger_df = fread(cellranger_metrics_file1)
    colnames(cellranger_df) = gsub(" ","_",colnames(cellranger_df))
    cellranger_ncells = as.numeric(gsub(",","",cellranger_df$Estimated_Number_of_Cells))
    df_to_plot = rbind(df_to_plot, data.frame(sample=sample, cellranger_ncells=cellranger_ncells))
}

df_to_plot1 = df_to_plot[1:120,]
df_to_plot2 = df_to_plot[121:nrow(df_to_plot),]

options(repr.plot.width = 18, repr.plot.height = 6)
p = ggplot(df_to_plot1, aes(x=sample, y=cellranger_ncells, fill=sample)) + geom_bar(stat = "identity")
p = p + theme_classic() + theme(text = element_text(size=12))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
p1 = p + ylab("Total cells") + ggtitle("CellRanger number of cells (TOB)")
p1

p = ggplot(df_to_plot2, aes(x=sample, y=cellranger_ncells, fill=sample)) + geom_bar(stat = "identity")
p = p + theme_classic() + theme(text = element_text(size=12))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
p2 = p + ylab("Total cells") + ggtitle("CellRanger number of cells (BioHEART)")
p2

# save CellRanger numbers
df_to_save = df_to_plot
df_to_save$cellbender_ncells <- c()
df_to_save$cohort = "TOB"
df_to_save$cohort[121:nrow(df_to_save)] = "BioHEART"
fwrite(df_to_save, "/share/ScratchGeneral/anncuo/tenk10k/data_processing/cellranger_ncells_summary.csv")

# save map from seq lib to cohort
df1 = df
colnames(df1)[1] = "sequencing_library"
df1$cellranger_ncells <- c()
fwrite(df1, "/share/ScratchGeneral/anncuo/tenk10k/data_processing/sequencing_library_to_cohort_map.csv")
