# Demultiplexing and doublet detection strategy

We use demuxafy ([preprint](https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1), [docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html)) to determine the tools to use for doublet detection as well as demultiplexing cells back to their donor of origin.

The scripts in this folder are structured as follows:

* Demultiplexing tools:
  * [Demuxalot script](demuxafy_demuxalot_runner.qsub) using [Demuxalot](https://pypi.org/project/demuxalot/), requires genotypes
  * [Dropulation script](demuxafy_dropulation_runner.qsub) using [Dropulation](https://github.com/broadinstitute/Drop-seq/blob/master/doc/Census-seq_Computational_Protcools.pdf), requires genotypes
  * [Vireo script](demuxafy_vireo_runner.qsub) using [vireo](https://vireosnp.readthedocs.io/en/latest/index.html)
* Doublet detection tools:
  * [Scds script](demuxafy_scds_runner.qsub) using [scds](https://github.com/kostkalab/scds)
  * [ScDblFinder script](demuxafy_scdblfinder_runner.qsub) uses [scDblFinder](https://github.com/plger/scDblFinder)
  * [DoubletDetection script](demuxafy_doubletdetection_runner.qsub) uses [DoubletDetection](https://github.com/JonathanShor/DoubletDetection)
* [Combiner script](demuxafy_combiner.qsub) combines all results using majority voting
