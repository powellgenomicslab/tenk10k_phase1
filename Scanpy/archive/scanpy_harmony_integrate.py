import glob
import matplotlib.pyplot as plt

import scanpy as sc
import scanpy.external as sce

# Filtered object directory
scanpy_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/filtered_scanpy_objects/"

# load adata
scanpy_files = glob.glob(scanpy_dir+"S*")

datasets=[]
for file in scanpy_files:
    adata = sc.read(file)
    datasets.append(adata)

adata = datasets[0].concatenate(*datasets[1:])

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


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
plt.savefig(scanpy_dir+"figures/merged_umap_leiden.pdf")

sc.pl.umap(adata, color=['sample'])
plt.savefig(scanpy_dir+"figures/merged_umap_sample.pdf")

sc.tl.tsne(adata)
sc.pl.tsne(adata, color=['leiden'])
plt.savefig(scanpy_dir+"figures/merged_tsne_leiden.pdf")

sc.pl.tsne(adata, color=['sample'])
plt.savefig(scanpy_dir+"figures/merged_tsne_sample.pdf")

sce.pp.harmony_integrate(adata, 'sample')

adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata)
sc.tl.tsne(adata)
sc.tl.leiden(adata, resolution=0.5)

out_dir = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/integration/"
out_file = out_dir+"harmony_TOB_samples_concatenated.h5ad"
adata.write(out_file)

sc.pl.umap(adata, color=['leiden'])
plt.savefig(scanpy_dir+"figures/merged_umap_leiden_harmony.pdf")

sc.pl.tsne(adata, color=['leiden'])
plt.savefig(scanpy_dir+"figures/merged_tsne_leiden_harmony.pdf")

sc.pl.umap(adata, color=['sample'])
plt.savefig(scanpy_dir+"figures/merged_umap_sample_harmony.pdf")

sc.pl.tsne(adata, color=['sample'])
plt.savefig(scanpy_dir+"figures/merged_tsne_sample_harmony.pdf")




