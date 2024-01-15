args <- commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

library(Seurat)

# 64 samples from 231013
# cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"

# 24 samples from 231213
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231213_tenk10k_gencode44/cellranger_outs/"

# 41 samples from 231214
# cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231214_tenk10k_gencode44/cellranger_outs/"

# 18 samples from 240108
cellranger_dir = "/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240108_tenk10k_gencode44/cellranger_outs/"

cellranger_samples = list.files(cellranger_dir, pattern = "S.") 

out_folder = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/seurat_objects/"

sample = cellranger_samples[i]

print(paste0("Create Seurat object for ",sample))

# cellranger_path = paste0(cellranger_dir,sample,"/outs/filtered_feature_bc_matrix/")
cellranger_path = paste0(cellranger_dir,sample,"/cellranger_count/",sample,"/outs/filtered_feature_bc_matrix/")
pbmc.data <- Read10X(data.dir = cellranger_path)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = sample, min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
out_file = paste0(out_folder, sample, ".rds")
saveRDS(pbmc, file = out_file)
