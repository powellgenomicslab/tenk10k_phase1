# Input files recap

## Common inputs

All tools require count information which gets outputted directly from Cellranger:

* scds, scDblFinder require **counts matrix (filtered)**: ```outs/filtered_feature_bc_matrix.h5```
* vireo requires **bams**: ```outs/possorted_genome_bam.bam```
* vireo also requires to specify **barcodes** (cell list), which can be extracted directly from cellranger (```outs/filtered_feature_bc_matrix/barcodes.tsv.gz```), although we will use those from cellbender (```cellbender_output_cell_barcodes.csv```)

