import scanpy as sc

sc.settings.verbosity = 2  # show logging output
sc.settings.autosave = True  # save figures, do not show them
sc.settings.set_figure_params(dpi=300)  # set sufficiently high resolution for saving
sc.settings.figdir = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures'
# filename = sys.argv[1]  # read filename from command line
filename = '/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info.h5ad'

def basic_analysis(filename):
    adata = sc.read_10x_h5(filename)
    sc.pp.recipe_zheng17(adata) # run the zheng17 recipe - used to process 1.3 million cells
    sc.pp.neighbors(adata)
    sc.tl.louvain(adata)
    sc.tl.paga(adata)
    sc.tl.umap(adata)
    sc.tl.rank_genes_groups(adata, 'louvain')
    adata.write('/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries_concatenated_gene_info_umaps.h5ad')
    # plotting
    sc.pl.paga(adata, save='.png')
    sc.pl.umap(adata, color='louvain', save='.png')
    sc.pl.rank_genes_groups(adata, save='.png')


if __name__ == "__main__":
    basic_analysis(filename)