library(data.table)
library(ggplot2)

demuxafy_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/"

# Vireo run using cellranger BAMs, including the number of expected donors even if not in the VCF file
vireo_dir_normal = paste0(demuxafy_dir, "vireo_output_nogeno_donors/")

# Vireo run with similar settings, but using the barcode file from Cellbender
vireo_dir_cellbender_barcodes = paste0(demuxafy_dir, "vireo_output_cellbender_barcodes/")

# Extract samples
samples = list.files(vireo_dir_cellbender_barcodes)

# Summarise first set of results
df_summary = data.frame()
for (sample in samples){
    file = paste0(vireo_dir_normal,sample,"/summary.tsv")
    vireo_df = read.csv(file, sep="\t")
    new_df = data.frame(sample = sample,
                        n_donors = nrow(vireo_df),
                        n_cells = sum(vireo_df$Freq),
                        n_doublets = vireo_df[vireo_df$Var1=="doublet","Freq"],
                        n_unassigned = vireo_df[vireo_df$Var1=="unassigned","Freq"])
    df_summary = rbind(df_summary, new_df)
}
df_summary1 = df_summary

# and same for second set
df_summary = data.frame()
for (sample in samples){
    file = paste0(vireo_dir_cellbender_barcodes, sample, "/summary.tsv")
    if (file.exists(file) == FALSE){next}
    vireo_df = read.csv(file, sep="\t")
    new_df = data.frame(sample = sample,
                        n_donors = nrow(vireo_df),
                        n_cells = sum(vireo_df$Freq),
                        n_doublets = vireo_df[vireo_df$Var1=="doublet","Freq"],
                        n_unassigned = vireo_df[vireo_df$Var1=="unassigned","Freq"])
    df_summary = rbind(df_summary, new_df)
}
df_summary2 = df_summary

# combine results
df_summary1$sample_cells = paste0(df_summary1$sample,"_cr")
df_summary2$sample_cells = paste0(df_summary2$sample,"_cb")
df_summary1$cells = "cellranger"
df_summary2$cells = "cellbender"
df_summary = rbind(df_summary1[df_summary1$sample %in% df_summary2$sample, ], df_summary2)

# define "assigned" cells as those assigned to a donor, so neither doublets nor unassigned
df_summary$n_assigned = df_summary$n_cells-df_summary$n_doublets-df_summary$n_unassigned


# plot assigned cells by sample + method
p = ggplot(df_summary, aes(x=sample_cells, y=n_assigned, alpha=cells, col=sample, fill=sample)) 
p = p + geom_bar(stat = "identity")
p = p + theme_classic() + theme(text = element_text(size=20))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
p = p + ylab("Total assigned cells") + ggtitle("vireo assigned (neither doublets nor unassigned)")
p 


# doublet rate
p = ggplot(df_summary, aes(x=sample_cells, y=n_doublets/n_cells, alpha=cells, col=sample, fill=sample)) 
p = p + geom_bar(stat = "identity")
p = p + theme_classic() + theme(text = element_text(size=20))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")
p = p + ylab("Doublet rate") + ggtitle("vireo doublets rates")
p


# unassigned rate
p = ggplot(df_summary, aes(x=sample_cells, y=n_unassigned/n_cells, alpha=cells, col=sample, fill=sample)) 
p = p + geom_bar(stat = "identity")
p = p + theme_classic() + theme(text = element_text(size=20))
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")
p = p + ylab("Unassigned rate") + ggtitle("vireo unassigned rate")
p
