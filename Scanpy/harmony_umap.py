import scanpy as sc
import scanpy.external as sce
# import matplotlib
# matplotlib.use('Agg')  # plotting backend compatible with screen

sc.settings.verbosity = 4  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=400)  # set sufficiently high resolution for saving
sc.settings.figdir = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures'
# filename = sys.argv[1]  # read filename from command line
filename = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info_harmony_umap_notsne.h5ad'

def basic_analysis(filename):

    adata = sc.read(filename) 

    # only consider genes with more than 3 counts
    sc.pp.filter_genes(adata, min_counts=3)  
    
    # Use Anna's preprocessing steps from prepare_cell_covs.py
    # sc.pp.normalize_total(adata, target_sum=1e4)
    # sc.pp.log1p(adata)
    # sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # adata.raw = adata
    # adata = adata[:, adata.var.highly_variable]
    # sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    # sc.pp.scale(adata, max_value=10)

    # Dim Reduction - without Harmony batch correction for comparison
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata) 
    sc.tl.umap(adata) 

    # plotting (No Harmony)
    sc.pl.pca(adata, color = 'wg2_scpred_prediction', save='_cell_types_noharmony.png')
    sc.pl.pca(adata, color = ['cohort'], save='pca_cohort_noharmony.png')
    
    # UMAP plotting (No Harmony)
    sc.pl.umap(adata, color='wg2_scpred_prediction', save='_celltypes_noharmony.png')
    sc.pl.umap(adata, color='cohort', save='_cohort_noharmony.png')

    # perform Harmony integration to remove technical batch effects between pools
    # sce.pp.harmony_integrate(adata, 'sequencing_library')

    # # Re-run UMAP + TSNE on the Harmony principal components 
    # adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
    # sc.pp.neighbors(adata)
    # sc.tl.umap(adata)

    # # plotting (Harmony Batch-corrected)
    # sc.pl.pca(adata, color = 'wg2_scpred_prediction', save='pca_cell_types.png')
    # sc.pl.pca(adata, color = ['cohort'], save='pca_cohort.png')
    # sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='')
    
    # # UMAP plotting (No Harmony)
    # sc.pl.umap(adata, color='louvain', save='umap_clusters.png')
    # sc.pl.umap(adata, color='wg2_scpred_prediction', save='umap_celltypes.png')
    # sc.pl.umap(adata, color='cohort', save='umap_cohort.png')

    # adata.write('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info_harmony_umap.h5ad')


if __name__ == "__main__":
    basic_analysis(filename)