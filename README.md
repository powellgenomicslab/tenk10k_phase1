# TenK10K phase1

Repo for analyses of the pilot phase of TenK10K.

## Pre-processing

* [CellBender](https://cellbender.readthedocs.io/en/latest/tutorial/index.html) for ambient RNA detection
* [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html) for demultiplexing and doublet detection
* QC & normalisation using [Scanpy](https://scanpy.readthedocs.io/en/stable/)
* batch correction / integration testing different strategies, mainly [BBKNN](https://github.com/Teichlab/bbknn) and [Harmony](https://portals.broadinstitute.org/harmony/) atm
* cell typing testing different strategies, starting with [CellTypist](https://www.celltypist.org/)

## Main analysis

* [SAIGE-QTL](https://github.com/weizhou0/qtl) for association between genetic variants (common and rare) and single-cell expression
* variant annotation
