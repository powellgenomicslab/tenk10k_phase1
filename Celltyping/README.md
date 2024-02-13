# Cell typing, _i.e._, assign cells to specific types of cells

## Reference-based approaches

For well-characterised biological systems such as circulating immune cells, or peripheral blood mononuclear cells (PBMCs) as we have here, we can rely on reference-based approaches, that use a reference dataset (often protein-based) to classify each cell into a known (sub) cell type, based on its expression profile.
The advantage of these approaches is that this analysis can be done one cell at a time, and therefore be run one batch / sequencing library at a time as opposed to requiring a fully integrated (potentially very large) dataset.

Here, we run three reference-based approaches:

* (hierarchical) scPred as recommended by the single-cell eQTL Gen consortium, see below
* Azimuth, also implemented as part of the single-cell eQTL Gen consortium
* celltypist

### scPred & Azimuth

The single-cell eQTL Gen consortium recommends the following approach, [instruction docs](https://powellgenomicslab.github.io/WG2-pipeline-classification-docs/general.html) for working group 2 (WG2), classification.

As recommended, we run these two tools using the singularity image provided.
As the input files expected here are Seurat objects (one per sequencing library), the first step here includes creating Seurat objects (directly from the CellRanger outputs), taking care not to include any QC steps here, so that the input objects are comparable to other tools used here

### CellTypist



## Brief description of scripts

* Create Seurat objects, shell script to run R script for each sequencing library independently
* Run hierarchical scPred shell script
* Run Azimuth shell script
* celltypist
