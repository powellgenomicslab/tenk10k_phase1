# QC and processing using Scanpy

We use Python's [Scanpy](https://scanpy.readthedocs.io/en/stable/) as main processing tool, as it scales better than R's Seurat / SingleCellExperiment.

## Scripts

* First, we add Cellbender, Demuxafy and cell typing results to each sequencing library's AnnData object, using [this script](), ran in parallel by [this runner]()
* Next, we perform basic pre-processing and QC, also one sequencing library at a time. [Script](), [runner]().
* Finally, we merge all objects into a single one, perform batch correction using Harmony, and plot the data. [Script]().


## Resources

* [video tutorial on use of Scanpy for scRNA-seq data](https://www.youtube.com/watch?v=uvyG9yLuNSE), found by [Tess](https://www.katalog.uu.se/profile/?id=N18-736)
