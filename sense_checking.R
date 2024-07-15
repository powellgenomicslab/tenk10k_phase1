library(Seurat)
library(tidyverse)

S001A_unfiltered <- readRDS("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/seurat_objects/S0001a.rds")
metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_cell_metadata_subset.csv")


S001A_barcodes_modified <- colnames(S001A_unfiltered) %>% str_replace("-1", "_S0001a") # update the barcode formatting

metadata <- metadata %>%
    rename(barcode = `...1`)
metadata$barcode %>% head()

length(setdiff(S001A_barcodes_modified, metadata$barcode)) # 6553 barcodes that were removed (likely outside QC thresholds)
length(intersect(S001A_barcodes_modified, metadata$barcode)) # 24896 barcodes kept
# this indicates that there was no issues with barcodes getting mixed up or something like that...

removed_barcodes_medatata <- S001A_unfiltered@meta.data[setdiff(S001A_barcodes_modified, metadata$barcode) %>% str_replace("_S0001a", "-1"), ]

head(removed_barcodes_medatata)
