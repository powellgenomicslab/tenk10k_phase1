# TenK10K phase1 scRNA-seq processing steps

Repo for analyses of the pilot phase of TenK10K (scRNA-seq processing only for now).

## Pipeline

* Create a file summarising sequencing libraries included in the latest run (by running [this](https://github.com/powellgenomicslab/tenk10k_phase1/blob/main/Demuxafy/preprocessing/prepare_inputs.md#extract-cellranger-libraries))
* Ambient RNA
  * [CellBender runner](CellBender/cellbender_runner.qsub), qsub script to run CellBender for each sequencing library
* Doublet Detection + Demultiplexing  
  * [scds runner](Demuxafy/demuxafy_scds_runner.qsub), qsub script to run scds for each sequencing library using Demuxafy image
  * [scDblFinder runner](Demuxafy/demuxafy_scdblfinder_runner.qsub), qsub script to run scDblFinder for each sequencing library using Demuxafy image
  * [vireo runner](Demuxafy/demuxafy_vireo_runner.qsub), qsub script to run vireo for each sequencing library using Demuxafy image (requires CellBender results, genotype info)
  * [Demuxafy combiner script](Demuxafy/demuxafy_combiner.qsub), qsub script to run Demuxafy combiner for each sequencing library (requires scds, scDblFinder, vireo results)
* Cell Typing
  * Consortium WG2 (scPred + Azimuth)  
    * [make Seurat objects scripts](Celltyping/make_Seurat_objects.qsub), qsub script running the [R script](Celltyping/make_Seurat_object_one_sample.R) building Seurat objects for each sequencing library prior to cell typing using Azimuth / scPred
    * [Azimuth cell typing runner](Celltyping/WG2_map_azimuth.qsub), qsub script to perform celltyping using Azimuth for each sequencing library using sceQTLGen WG2 image (requires Seurat objects)
    * [hierarchical scPred cell typing runner](Celltyping/WG2_map_hierscpred.qsub), qsub script to perform celltyping using hierarchical scPred for each sequencing library using sceQTLGen WG2 image (requires Seurat objects)
    * [Consortium WG2 cell typing combiner script](Celltyping/combine_results.qsub), qsub script running the [R script](Celltyping/combine_results_one_sample.R) which combines cell types prediction for each sequencing library obtained using Azimuth and scPred (requires azimuth and hierarchical scPred results)
  * Celltypist
    * [Celltypist runner](Scanpy/celltypist_runner.qsub), qsub script running the [Python script](Scanpy/celltypist_per_sample.py) performing cell typing using Cell Typist for each sequencing library
* Scanpy data wrangling & data integration
  * [Add info runner](Scanpy/run_add_metadata.qsub), qsub script running the [Python script](Scanpy/add_metadata_per_sample.py) adding all metadata to scanpy object for each sequencing library (requires results from CellBender, Demuxafy combiner, WG2 cell typing combiner, celltypist results) and performing initial QC
  * Python script combining results into a single AnnData object and performing batch correction + data integration using Harmony 

## Pre-processing tools used

* [CellBender](https://cellbender.readthedocs.io/en/latest/tutorial/index.html) for ambient RNA detection
* [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html) for demultiplexing and doublet detection, specifically:
  * [vireo](https://vireosnp.readthedocs.io/en/latest/manual.html) for demultiplexing, as it allows to include the number of donors expected per pool, regardless of whether we have genotype data for them
  * majority voting of [vireo](https://vireosnp.readthedocs.io/en/latest/manual.html), [scds](https://github.com/kostkalab/scds) and [scDblFinder](https://github.com/plger/scDblFinder) for doublet detection 
* QC & normalisation using [Scanpy](https://scanpy.readthedocs.io/en/stable/)
* batch correction / integration using [Harmony](https://portals.broadinstitute.org/harmony/) 
* cell typing using [scPred](), [Azimuth](), and [CellTypist](https://www.celltypist.org/)

## Results on Brenner

* CellBender results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/cellbender_output_smaller_learning_rate/{SAMPLE}/```
* scds results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/scds_output/{SAMPLE}/```
* scDblFinder results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/scdblfinder_output/{SAMPLE}/```
* vireo results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/vireo_output_final/{SAMPLE}/```
* Consortium WG2 Azimuth results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/{SAMPLE}/step2_azimuth/```
* Consortium WG2 hierarchical scPred results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/{SAMPLE}/step3_hierscpred/```

## Keeping track of samples generated (scRNA-seq)

* 1st batch: Oct 13th, 2023 (aka ```231013```): 64 sequencing libraries (TOB only), 578 unique individuals
* 2nd batch: Dec 13th, 2023 (```231213```): 24 sequencing libraries (TOB only), 222 unique individuals
* 3rd batch: Dec 14th, 2023 (```231214```): 41 sequencing libraries (24 TOB, 17 BioHEART), 343 unique individuals (210 TOB, 133 BioHEART).
* 4th batch: Jan 8th, 2024 (```240108```): 18 sequencing libraries (BioHEART only), 140 unique individuals.
* 5th batch: Jan 12th, 2024 (```240112```): 18 sequencing libraries (BioHEART only), 139 unique individuals.
* 6th batch: Jan 15th, 2024 (```240115```): 25 sequencing libraries (8 TOB, 17 BioHEART), 173 unique individuals (112 TOB, 61 BioHEART).
* 7th batch: Jan 16th, 2024 (```240116```): 17 sequencing libraries (BioHEART only), 132 unique individuals.
* 8th batch: Jan 19th, 2024 (```240119```): 17 sequencing libraries (BioHEART only), XX unique individuals.
