import scanpy as sc
import scanpy.external as sce
# import matplotlib
# matplotlib.use('Agg')  # plotting backend compatible with screen

sc.settings.verbosity = 4  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi_save=500)  # set sufficiently high resolution for saving and so the colors dont get mixed wierdly
sc.settings.figdir = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures'
# filename = sys.argv[1]  # read filename from command line
filename = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info_harmony_umap_notsne.h5ad'

def basic_analysis(filename):

    adata = sc.read(filename) 
    # Re-run UMAP + TSNE on the Harmony principal components
    # adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
    # sc.pp.neighbors(adata)
    # sc.tl.umap(adata)
    # sc.tl.tsne(adata)

    # PCA plotting (Harmony Batch-corrected)
    sc.pl.pca(adata, color = 'wg2_scpred_prediction', save='pca_cell_types2.png')
    sc.pl.pca(adata, color = ['cohort'], save='pca_cohort2.png')
    
    # UMAP plotting Harmony batch-corrected 
    sc.pl.umap(adata, color='wg2_scpred_prediction', save='umap_celltypes2.png')
    sc.pl.umap(adata, color='cohort', save='umap_cohort2.png')

    # adata.write('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info_harmony_umap_notsne.h5ad')

if __name__ == "__main__":
    basic_analysis(filename)