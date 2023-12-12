import glob
import scanpy as sc
import scanpy.external as sce

import celltypist
from celltypist import models

# Filtered object directory
scanpy_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/filtered_scanpy_objects/"

# load adata
scanpy_files = glob.glob(scanpy_dir+"S*")

datasets=[]
for file in scanpy_files:
    adata = sc.read(file)
    datasets.append(adata)

adata = datasets[0].concatenate(*datasets[1:])

model = models.Model.load(model = 'Immune_All_Low.pkl')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_pcs = 30)

sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 0.5)
sc.pl.umap(adata, color=['leiden'])

sc.pl.umap(adata, color=['sample'])

sc.tl.tsne(adata)
sc.pl.tsne(adata, color=['leiden'])

sc.pl.tsne(adata, color=['sample'])

sce.pp.harmony_integrate(adata, 'sample')

adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.tsne(adata)
sc.tl.leiden(adata, resolution=0.5)

sc.pl.umap(adata, color=['leiden'])
sc.pl.tsne(adata, color=['leiden'])

sc.pl.umap(adata, color=['sample'])
sc.pl.tsne(adata, color=['sample'])

out_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/integration/"
out_file = out_dir+"harmony_64_samples_concatenated.h5ad"
adata.write(out_file)


