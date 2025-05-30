import pandas as pd
import scanpy as sc

# AnnData objects location, split by chrom + cell type
anndata_dir = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/cpg_anndata/'

# output directory
out_dir = '/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/review_files/ct_mean_expression_files/'

# extract unique cell type from eQTL summary file
ct_df = pd.read_csv('/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/ms_tables/all_eqtls_fdr_5pct.tsv')
celltypes = ct_df.celltype.unique()

# loop over cell types and chromosomes to open each AnnData object
for celltype in celltypes:
  for chrom in range(22):
    chrom=chrom+1
    # open file
    adata_file = f'{anndata_dir}{celltype}_chr{chrom}.h5ad'
    adata = sc.read(adata_file)
    # extract mean expression across all cells for a given cell type (wg2_scpred_prediction)
    means_by_ct = sc.get.aggregate(adata, by="wg2_scpred_prediction", func="mean")
    # save as pandas data frame
    expr_df = pd.DataFrame(data=means_by_ct.layers['mean'], index=[f'mean_{celltype}'], columns = means_by_ct.var.index)
    # write to file
    print('Writing to file: {celltype}, {chrom}')
    out_file = f'{out_dir}{celltype}_chr{chrom}.csv'
    expr_df.to_csv(out_file)
