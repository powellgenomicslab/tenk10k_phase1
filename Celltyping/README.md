# Cell typing, i.e. assign cells to specific types of cells

## Reference-based approaches

For well-characterised biological systems such as circulating immune cells, or peripheral blood mononuclear cells (PBMCs) as we have here, we can rely on reference-based approaches, that use a reference dataset (often protein-based) to classify each cell into a known (sub) cell type, based on its expression profile.
The advantage of these approaches is that this analysis can be done one cell at a time, and therefore be run one batch / sequencing library at a time as opposed to requiring a fully integrated (potentially very large) dataset.

Reference-based approaches include:

* (hierarchical) scPred
* celltypist

### scPred

The single-cell eQTL Gen consortium recommends the following approach, [instruction docs](https://powellgenomicslab.github.io/WG2-pipeline-classification-docs/general.html) for working group 2 (WG2), classification.

* download singularity image
* create Seurat object (one per sequencing library)
* run map_azimuth.R
* run map_hierscpred.R
* compare
* maybe merge?
