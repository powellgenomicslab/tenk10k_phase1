# Demultiplexing and doublet detection strategy

We use demuxafy ([preprint](https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1), [docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html)) to determine the tools to use for doublet detection as well as demultiplexing cells back to their donor of origin.

The scripts in this folder are structured as follows:

* [Demuxalot script]()
* [Dropulation script]()
* [Scds script](demuxafy_scds_runner.qsub)
* [ScDblFinder script]()
* [DoubletDetection script]()
* [Combiner script](demuxafy_combiner.qsub) using majority voting
