# Input files recap

## Common inputs (count data from CellRanger)

All tools require count information which gets outputted directly from Cellranger:

* scds, scDblFinder require **counts matrix (filtered)**: ```outs/filtered_feature_bc_matrix.h5```
* vireo requires **bams**: ```outs/possorted_genome_bam.bam```
* vireo also requires to specify **barcodes** (cell list), which can be extracted directly from cellranger (```outs/filtered_feature_bc_matrix/barcodes.tsv.gz```), although we will use those from cellbender (```cellbender_output_cell_barcodes.csv```)

## Extra inputs for vireo (genotype and individual info)

Additionally, in order to demultiplex (assign cells to donors, _e.g._, using vireo), genotype files are required.

Vireo requires:

* a **VCF** file including genotypes for all individuals (this needs to be zipped and indexed, as well as other possible manipulations, see [this paragraph in the preprocessing](https://github.com/powellgenomicslab/tenk10k_phase1/blob/main/Demuxafy/preprocessing/prepare_inputs.md#vcf-manipulation))
* to specify the **individuals** from the VCF that are expected in the specific pool / library analysed
  * use all of the individuals [extracted from the VCF](https://github.com/powellgenomicslab/tenk10k_phase1/blob/main/Demuxafy/preprocessing/prepare_inputs.md#individual-files-directly-from-vcf)
  * or specify the specific ones included in the pool e.g. using [this script](TOB_make_pool_sample_lists.R) 
* it is possible to specify the total **number of individuals** in the pool (even if not in the VCF) using the -N flag, for which we [created number of individual files](TOB_make_pool_number_of_sample_lists.R)
