library(Seurat)

# 64 samples from 231013
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
cellranger_samples = list.files(cellranger_dir, pattern = "S.") 

out_folder = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/scpred/seurat_objects/"

for (sample in cellranger_samples){
  cellranger_path = paste0(cellranger_dir,sample,"/outs/filtered_feature_bc_matrix/")
  pbmc.data <- Read10X(data.dir = cellranger_path)
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = sample, min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  out_file = paste0(out_folder, sample, ".rds")
  saveRDS(pbmc, file = out_file)
}


