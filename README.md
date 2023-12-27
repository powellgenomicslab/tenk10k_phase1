# TenK10K phase1

Repo for analyses of the pilot phase of TenK10K.

## Pipeline

* [CellBender runner](CellBender/cellbender_runner.qsub), qsub script to run CellBender for each sequencing library 
* [scds runner](Demuxafy/demuxafy_scds_runner.qsub), qsub script to run scds for each sequencing library using Demuxafy image
* [scDblFinder runner](Demuxafy/demuxafy_scdblfiner_runner.qsub), qsub script to run scDblFinder for each sequencing library using Demuxafy image
* [vireo runner](Demuxafy/demuxafy_vireo_runner.qsub), qsub script to run vireo for each sequencing library using Demuxafy image (requires CellBender results, genotype info)
* [Demuxafy combiner script](Demuxafy/demuxafy_combiner.qsub), qsub script to run Demuxafy combiner for each sequencing library (requires scds, scDblFinder, vireo results)
* 

## Pre-processing

* [CellBender](https://cellbender.readthedocs.io/en/latest/tutorial/index.html) for ambient RNA detection
* [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html) for demultiplexing and doublet detection, specifically:
  * [vireo] for demultiplexing, as it allows to include the number of donors expected per pool, regardless of whether we have genotype data for them
  * majority voting of [vireo], [scds] and [scDblFinder] for doublet detection 
* QC & normalisation using [Scanpy](https://scanpy.readthedocs.io/en/stable/)
* batch correction / integration using [Harmony](https://portals.broadinstitute.org/harmony/) 
* cell typing testing different strategies, starting with [CellTypist](https://www.celltypist.org/)

## Main analysis

* [SAIGE-QTL](https://github.com/weizhou0/qtl) for association between genetic variants (common and rare) and single-cell expression
* variant annotation
