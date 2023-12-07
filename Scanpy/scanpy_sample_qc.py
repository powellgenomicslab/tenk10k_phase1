import glob
import sys
import scanpy as sc

i = int(sys.argv[1])

# CellRanger files
cellranger_dir = "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/"
cellranger_files = glob.glob(cellranger_dir+"S*")

sample_path = cellranger_files[i]

# load adata
adata=sc.read_10x_h5(sample_path+"/outs/filtered_feature_bc_matrix.h5")

print(f'Total number of cells: {adata.shape[0]}')
print(f'Total number of genes: {adata.shape[1]}')

# unique gene names
adata.var_names_make_unique()

print(f'Total number of unique genes: {adata.shape[1]}')

# basic filtering of cells and genes (should this be done when merged?)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# define MT genes and calculate QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# plot
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# apply filters based on plots
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]

# normalise
sc.pp.normalize_total(adata, target_sum=1e4) #normalize every cell to 10,000 UMI
sc.pp.log1p(adata)

# calculate and plot highly variable genes (HVGs)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)

# save all data as raw and proceed with HVGs only
adata.raw = adata
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'pct_counts_ribo'])
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_pcs = 30)
