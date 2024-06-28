# Plotting notebooks

Notebooks to reproduce figures generated in the analysis of these data

## [QC](qc_plots)

These are plots summarising the results of various scRNA-seq processing tools applied on this cohort, including ambient RNA and doublet detection, demultiplexing etc.

Most of these plots are split by sequencing library as most of these tools were run on each separately.

The [QC plots README](qc_plots/README.md) contains more detail. 

## [Data overview](overview_figures)

Thes are plots summarising the scRNA-seq data overall (across all sequencing libraries) and will constitute panels for Figure 1 and accompanying supplementary figures.

The [scRNA-seq data overview README](overview_figures/README.md) contains more detail.

### Main

* UMAP plot of all data coloured by cell type
* barplot of number of cells per cell type coloured by cell type
* barplot of number of cells per cell type coloured by cell type (log 10 scale)
* violin plots of number of cells per individual by cell type coloured by cell type
* histogram of number of cells per individual (across all cell types)

### Supplement

* QC metrics by cell type (violin plots coloured and grouped by cell type)
  * total counts
  * % MT 
