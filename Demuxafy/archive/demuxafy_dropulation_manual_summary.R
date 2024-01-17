#!/usr/bin/env Rscript
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]

fdrThreshold = 0.05
doubletPvalueThreshold = 1

results_dir = paste0("/share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/dropulation_output/",sample)

assign = paste0(results_dir, "/assignments.tsv.gz")
doubletFile = paste0(results_dir,"/likelihoods.tsv.gz")
out = paste0(results_dir,"/manual_summary.tsv")

a <- read.csv(gzfile(assign), sep="\t",header = FALSE)
colnames(a) = a[2,]
a = a[3:nrow(a),]
a$FDR_pvalue = as.numeric(a$FDR_pvalue)

b=fread(doubletFile, stringsAsFactors = F, sep="\t")

confidentAssignmentCells=a[a$FDR_pvalue<=fdrThreshold,"cell"]
singletCells=b[b$doublet_pval<=doubletPvalueThreshold, "cell"]$cell
confidentAssignmentSingletCells=intersect(confidentAssignmentCells, singletCells)
doublets = a$cell[!(a$cell %in% confidentAssignmentSingletCells)]

mapa=a[match(confidentAssignmentSingletCells, a$cell), c("cell", "bestLikelihood", "bestSample")]
colnames(mapa) = c("Barcode","dropulation_Likelihood", "dropulation_Assignment")
mapa$dropulation_DropletType = "singlet"

mapb= b[match(doublets, b$cell), c("cell", "mixedSampleLikelihood", "mixedSample")]
colnames(mapb) = c("Barcode", "dropulation_Likelihood","dropulation_Assignment")
mapb$dropulation_Assignment = "doublet"
mapb$dropulation_DropletType = "doublet"

map = rbind(mapa,mapb)

# left join
a = data.table(a)
a$Barcode = a$cell
map = a[,c("Barcode", "num_snps", "num_umis")][map, on = "Barcode"]

colnames(map) = gsub("num_snps", "dropulation_Nsnps", colnames(map)) %>%
                        gsub("num_umis", "dropulation_Numis", .)

fwrite(map, out, row.names=F, col.names = T, quote=F, sep="\t")
