import glob
import scanpy as sc

# CellRanger files
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
files = glob.glob(cellranger_dir+"S*")

# load adata
adata=sc.read_10x_h5(myfile+"/outs/filtered_feature_bc_matrix.h5")

datasets = []

for file in files:
    sample = file.replace(cellranger_dir, "")
    # adata = sc.read(file)
    adata=sc.read_10x_h5(file+"/outs/filtered_feature_bc_matrix.h5")
    adata.obs['sample'] = sample
    adata.obs['barcode'] = adata.index
    adata.index = sample + "_" + adata.obs['barcode']
    datasets.append(adata)

adata_all = datasets[0].concatenate(*datasets[1:])

output_dir = 
filename_all = output_dir+"all_samples_concatenated.h5"
adata_all.write(filename_all)
