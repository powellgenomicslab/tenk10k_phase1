library(data.table)
library(ggplot2)

summary_filename = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/cellbender/cells_set_celltypist_confidence_summary_uncorrected.csv"
df = read.csv(summary_filename, row.names = 1)
samples = df$samples

options(repr.plot.width = 15, repr.plot.height = 8) 
p = ggplot(df[df$samples %in% samples,], aes(x=samples, y=cellranger_avg_celltypist_confidence)) + geom_point(size=3) 
p = p + geom_point(aes(x=samples, y=cellbender_09_avg_celltypist_confidence), data = df[df$samples %in% samples,], col="blue", size=3)
p = p + geom_point(aes(x=samples, y=cellbender_05_avg_celltypist_confidence), data = df[df$samples %in% samples,], col="red", size=3)
p = p + theme_minimal() + theme(text = element_text(size=20)) + ylab("celltypist average confidence score")
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1 = p + ggtitle("uncorrected counts")


summary_filename = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/cellbender/cells_set_celltypist_confidence_summary_corrected.csv"
df = read.csv(summary_filename, row.names = 1)
samples = df$samples

options(repr.plot.width = 15, repr.plot.height = 8) 
p = ggplot(df[df$samples %in% samples,], aes(x=samples, y=cellranger_avg_celltypist_confidence)) + geom_point(size=3) 
p = p + geom_point(aes(x=samples, y=cellbender_09_avg_celltypist_confidence), data = df[df$samples %in% samples,], col="blue", size=3)
p = p + geom_point(aes(x=samples, y=cellbender_05_avg_celltypist_confidence), data = df[df$samples %in% samples,], col="red", size=3)
p = p + theme_minimal() + theme(text = element_text(size=20)) + ylab("celltypist average confidence score")
p = p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2 = p + ggtitle("corrected counts")

