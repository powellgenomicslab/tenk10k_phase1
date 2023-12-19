# Input files recap

## Common inputs (CellRanger)

All tools require count information which gets outputted directly from Cellranger:

* scds, scDblFinder require **counts matrix (filtered)**: ```outs/filtered_feature_bc_matrix.h5```
* vireo requires **bams**: ```outs/possorted_genome_bam.bam```
* vireo also requires to specify **barcodes** (cell list), which can be extracted directly from cellranger (```outs/filtered_feature_bc_matrix/barcodes.tsv.gz```), although we will use those from cellbender (```cellbender_output_cell_barcodes.csv```)

## Genotypes

Additionally, in order to demultiplex (assign cells to donors, _e.g._, using vireo), genotype files are required.

Vireo requires:

* a VCF file including genotypes for all individuals (this needs to be zipped and indexed, as well as other possible manipulations, see prepare_inputs.md#vcf-manipulation)
* to specify the individuals from the VCF that are expected in the specific pool / library analysed
* it is possible to specify the total number of individuals in the pool (even if not in the VCF) using the -N flag, for which we created number of individual files
