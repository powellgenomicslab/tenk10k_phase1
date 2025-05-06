![singlecell pipeline data flow (for sharing)](https://github.com/user-attachments/assets/298491c4-273e-49e6-baab-5ce0e61c4c92)# Impact of rare and common genetic variation on cell type-specific gene expression

This repository contains scripts for data processing, analysis and figure generation using data from phase 1 of the TenK10K project for our paper:

Cuomo et al., [Impact of rare and common genetic variation on cell type-specific gene expression](https://www.medrxiv.org/content/10.1101/2025.03.20.25324352v1), medRxiv, 2025.

This includes:

* scRNA-seq processing 
* [cell-state abundance QTLs](cell_state_abundance_qtl)
* [manuscript main and suppl figures](plotting_notebooks)

## scRNA-seq processing pipeline

* Create a file summarising sequencing libraries included in the latest run (by running [this](https://github.com/powellgenomicslab/tenk10k_phase1/blob/main/Demuxafy/preprocessing/prepare_inputs.md#extract-cellranger-libraries))
* New pool-donor info (for vireo) may not be in the files used [here](https://github.com/powellgenomicslab/tenk10k_phase1/blob/main/Demuxafy/preprocessing/BioHEART_make_pool_sample_lists.R#L3-L6), but they will be in the [shared data tracking spreadsheet](https://docs.google.com/spreadsheets/d/1XcAYd3ssKTfdPquQ3QxtN8Mn18YObpiiuUecQFbz45Q/edit#gid=336263672), so download the updated list from there.
* Ambient RNA
  * [CellBender runner](CellBender/cellbender_runner.qsub), qsub script to run CellBender for each sequencing library
* Doublet Detection + Demultiplexing  
  * [scds runner](Demuxafy/demuxafy_scds_runner.qsub), qsub script to run scds for each sequencing library using Demuxafy image
  * [scDblFinder runner](Demuxafy/demuxafy_scdblfinder_runner.qsub), qsub script to run scDblFinder for each sequencing library using Demuxafy image
  * [vireo runner](Demuxafy/demuxafy_vireo_runner.qsub), qsub script to run vireo for each sequencing library using Demuxafy image (requires CellBender results, genotype info) -- note that running vireo requires knowing which individuals we expect in each pool, detailed in scripts in [this folder](/Demuxafy/preprocessing/)
  * [Demuxafy combiner script](Demuxafy/demuxafy_combiner.qsub), qsub script to run Demuxafy combiner for each sequencing library (requires scds, scDblFinder, vireo results)
* Cell Typing
  * Consortium WG2 (scPred + Azimuth)  
    * [make Seurat objects scripts](Celltyping/scpred/make_Seurat_objects.qsub), qsub script running the [R script](Celltyping/scpred/make_Seurat_object_one_sample.R) building Seurat objects for each sequencing library prior to cell typing using Azimuth / scPred
    * [Azimuth cell typing runner](Celltyping/scpred/WG2_map_azimuth.qsub), qsub script to perform celltyping using Azimuth for each sequencing library using sceQTLGen WG2 image (requires Seurat objects)
    * [hierarchical scPred cell typing runner](Celltyping/scpred/WG2_map_hierscpred.qsub), qsub script to perform celltyping using hierarchical scPred for each sequencing library using sceQTLGen WG2 image (requires Seurat objects)
    * [Consortium WG2 cell typing combiner script](Celltyping/scpred/combine_results.qsub), qsub script running the [R script](Celltyping/combine_results_one_sample.R) which combines cell types prediction for each sequencing library obtained using Azimuth and scPred (requires azimuth and hierarchical scPred results)
  * Celltypist
    * [Celltypist runner](Celltyping/celltypist/run_celltypist.qsub), qsub script running the [Python script](Celltyping/celltypist/celltypist_per_sample.py) performing cell typing using Cell Typist for each sequencing library
* Scanpy data wrangling & data integration (to do AFTER everything else)
  * [Add info runner](Scanpy/run_add_metadata.qsub), qsub script running the [Python script](Scanpy/add_metadata_per_sample_no_norm.py) adding all metadata to scanpy object for each sequencing library (requires results from CellBender, Demuxafy combiner, WG2 cell typing combiner, celltypist results) and performing initial QC
  * Python script combining results into a single AnnData object (concatenate), and add gene and donor info ([Python script](Scanpy/combine_files_add_gene_info.py))
  * Python script making AnnData objects for each cell type + chromosome combination. [Python script](Scanpy/prepare_pheno.py), [qsub runner](Scanpy/prepare_pheno_runner.qsub).
  * Python script making TSV objects containing expression PCs for each cell type. [Python script](Scanpy/prepare_cell_covs.py), [qsub runner](Scanpy/prepare_covs_runner.qsub).
  * TO DO: Python script performing QC & processing on combined AnnData + plotting.

## Pre-processing tools used

* [CellBender](https://cellbender.readthedocs.io/en/latest/tutorial/index.html) for ambient RNA detection
* [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html) for demultiplexing and doublet detection, specifically:
  * [vireo](https://vireosnp.readthedocs.io/en/latest/manual.html) for demultiplexing, as it allows to include the number of donors expected per pool, regardless of whether we have genotype data for them
  * majority voting of [vireo](https://vireosnp.readthedocs.io/en/latest/manual.html), [scds](https://github.com/kostkalab/scds) and [scDblFinder](https://github.com/plger/scDblFinder) for doublet detection 
* QC & normalisation using [Scanpy](https://scanpy.readthedocs.io/en/stable/)
* batch correction / integration using [Harmony](https://portals.broadinstitute.org/harmony/) 
* cell typing using [scPred](https://github.com/powellgenomicslab/scPred), [Azimuth](https://satijalab.github.io/azimuth/index.html), and [CellTypist](https://www.celltypist.org/)

## Processing pipeline data flow diagram

![singlecell pipeline data flow (for sharing)](https://github.com/user-attachments/assets/5c415a2c-2931-4ec6-8750-3eda7293ddb7)



