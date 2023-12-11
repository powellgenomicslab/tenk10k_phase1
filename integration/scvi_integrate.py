import glob
import time
import numpy as np
import pandas as pd
import scanpy as sc
import scvi

# CellRanger files
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
files = glob.glob(cellranger_dir+"S*")

# pre-processing steps
# modify parameters based on sample-specific plots
def pp(sample_path):
    adata=sc.read_10x_h5(sample_path+"/outs/filtered_feature_bc_matrix.h5")
    sc.pp.filter_genes(adata, min_cells = 10)
#     sc.pp.highly_variable_genes(adata, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')
#     sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
     
    sample = file.replace(cellranger_dir, "")
    adata.var_names_make_unique()
    adata.obs['sample'] = sample
    adata.obs['barcode'] = adata.obs.index
    adata.obs.index = sample + "_" + adata.obs['barcode']
        
    sc.pp.filter_cells(adata, min_genes=200) #get rid of cells with fewer than 200 genes
    #sc.pp.filter_genes(adata, min_cells=3) #get rid of genes that are found in fewer than 3 cells
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim]
    adata = adata[adata.obs.pct_counts_mt < 15]

    return adata

datasets = []
for file in files:
    datasets.append(pp(file))
    
adata = datasets[0].concatenate(*datasets[1:])
sc.tl.pca(adata, svd_solver='arpack')
scvi.model.SCVI.setup_anndata(adata, layer = "counts",
                              categorical_covariate_keys=["sample"],
                              continuous_covariate_keys=['pct_counts_mt', 'total_counts'])
model = scvi.model.SCVI(adata)
model.train()
adata.obsm['X_scVI'] = model.get_latent_representation()
adata.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)
sc.pp.neighbors(adata, use_rep = 'X_scVI')
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.5)
sc.pl.umap(adata, color = ['leiden', 'Sample'], frameon = False)
adata.uns['scvi_markers'] = markers_scvi
adata.uns['markers'] = markers
adata.write_h5ad('integrated.h5ad')
model.save('model.model')

