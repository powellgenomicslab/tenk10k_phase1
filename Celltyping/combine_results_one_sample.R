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
if (file.exists(out_file) == TRUE){quit(save = "no", status = 1, runLast = FALSE)}

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
