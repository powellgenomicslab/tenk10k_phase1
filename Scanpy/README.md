# QC and processing using Scanpy

We use Python's [Scanpy](https://scanpy.readthedocs.io/en/stable/) as main processing tool, as it scales better than R's Seurat / SingleCellExperiment.

## Scripts

* First, we add Cellbender, Demuxafy and cell typing results to each sequencing library's AnnData object, using [this script](), run in parallel using [this runner]()

Next, prepare input files for association analysis (which we will transfer to GCP for SAIGE-QTL running):

* Concatenate all raw objects with meta, add gene info
* Split by chromosome, cell type -- this files contain most necessary information: gene expression (all genes), individual id and batch information, gene information to determine cis window

Separately, to build cell type expression covariates:

* One sequencing library at a time, we perform basic pre-processing and QC. [Script](), [runner]().
* We merge all objects into a single one, then split by covariates, integrate using Harmony, generate (Harmony) principal components, one cell type at a time

Finally, for visualisation purposes, 

* we merge all objects into a single one, perform batch correction using Harmony (all cell types), and plot the data. [Script]().


## Resources

* [video tutorial on use of Scanpy for scRNA-seq data](https://www.youtube.com/watch?v=uvyG9yLuNSE), found by [Tess](https://www.katalog.uu.se/profile/?id=N18-736)
