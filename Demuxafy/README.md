# Demultiplexing and doublet detection strategy

We use demuxafy ([preprint](https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1), [docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html)) by [Drew Neavin](https://github.com/drneavin) to determine the tools to use for doublet detection as well as demultiplexing cells back to their donor of origin.

## Description of scripts

The scripts in this folder are structured as follows:

### Runners

* Demultiplexing tools:
  * [Vireo script](demuxafy_vireo_runner.qsub) using [vireo](https://vireosnp.readthedocs.io/en/latest/index.html)
* Doublet detection tools:
  * [Scds script](demuxafy_scds_runner.qsub) using [scds](https://github.com/kostkalab/scds)
  * [ScDblFinder script](demuxafy_scdblfinder_runner.qsub) uses [scDblFinder](https://github.com/plger/scDblFinder)
* [Combiner script](demuxafy_combiner.qsub) combines all results using majority voting

### Pre-processing

* [preprocessing scripts folder](preprocessing/)

<!---
### Post-processing

* [demuxafy_dropulation_manual_summary.R](demuxafy_dropulation_manual_summary.R)

### other tools
* [Demuxalot script](demuxafy_demuxalot_runner.qsub) using [Demuxalot](https://pypi.org/project/demuxalot/), requires genotypes
* [Dropulation script](demuxafy_dropulation_runner.qsub) using [Dropulation](https://github.com/broadinstitute/Drop-seq/blob/master/doc/Census-seq_Computational_Protcools.pdf), requires genotypes
* [DoubletDetection script](demuxafy_doubletdetection_runner.qsub) uses [DoubletDetection](https://github.com/JonathanShor/DoubletDetection)
--->
