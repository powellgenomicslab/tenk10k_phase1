import scanpy as sc
import pandas as pd
import scanpy.external as sce

sc.settings.verbosity = 4  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=400, dpi_save=500)  # set sufficiently high resolution for saving
sc.settings.figdir = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures'
# filename = sys.argv[1]  # read filename from command line
filename = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony_leiden_filtered.h5ad'

def analysis_harmony_secondpass(filename):
    """ perform normalisation, scaling, feature selection, harmony integration to remove batch effects, and dimensionality reduction on the input scanpy object
        also generates a few UMAP and PCA figures 
    """

    adata = sc.read(filename) 
    adata = adata.raw.to_adata()

    # only consider genes with more than 3 counts across all cells 
    sc.pp.filter_genes(adata, min_counts=3)  

    # Use Anna's preprocessing steps from prepare_cell_covs.py 
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    # Dim Reduction - without Harmony batch correction for comparison
    sc.tl.pca(adata, svd_solver='arpack')

    # perform Harmony integration to remove technical batch effects between pools
    sce.pp.harmony_integrate(adata, 'sequencing_library')

    # Re-run UMAP on the Harmony principal components
    adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Save the object and metadata 
    adata.write('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_concatenated_harmony_leiden_filtered_reanalysed.h5ad')

    # save the UMAP coords for plotting in R
    pd.DataFrame(adata.obsm['X_umap'],
                columns = ['UMAP1', 'UMAP2'],
                index = adata.obs.index).to_csv('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_harmony_umap_coords_filtered_reanalysed.csv')

    # get a subset of the metadata for UMAP QC plots
    adata.obs[['wg2_scpred_prediction', 'leiden', 'individual', 'cohort', 'n_genes', 'n_genes_by_counts',
        'total_counts', 'total_counts_mt', 'pct_counts_mt', 'ct_id', 'cpg_id', 'pct_counts_ribo', 'pct_counts_hb']].to_csv('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_cell_metadata_subset_filtered_reanalysed.csv')


if __name__ == "__main__":
    analysis_harmony_secondpass(filename)


#  adata.obs[['wg2_scpred_prediction', 'cohort', 'individual']].to_csv('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_archive/240_libraries_celltypes_by_cohort_indiv.csv')
