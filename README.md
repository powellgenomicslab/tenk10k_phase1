# TenK10K phase1

Repo for analyses of the pilot phase of TenK10K.

## Pre-processing

* [CellBender](https://cellbender.readthedocs.io/en/latest/tutorial/index.html) for ambient RNA detection
* [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html) for demultiplexing and doublet detection
* QC & normalisation using [Scanpy](https://scanpy.readthedocs.io/en/stable/)
* batch correction / integration testing different strategies, mainly BBKNN and Harmony atm
* cell typing testing different strategies, starting with [CellTypist](https://www.celltypist.org/)

## Main analysis

* variant annotation
* SAIGE-QTL
