# Demultiplexing and doublet detection strategy

We use demuxafy ([preprint](https://www.biorxiv.org/content/10.1101/2022.03.07.483367v1), [docs](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html)) by [Drew Neavin](https://github.com/drneavin) to determine the tools to use for doublet detection as well as demultiplexing cells back to their donor of origin.

## Description of scripts

The scripts in this folder are structured as follows:

### Runners

* Demultiplexing tools:
  * [Vireo script (TOB)](demuxafy_vireo_runner_tob.qsub) using [vireo](https://vireosnp.readthedocs.io/en/latest/index.html)
  * [Vireo script (BioHEART)](demuxafy_vireo_runner_bioheart.qsub) using [vireo](https://vireosnp.readthedocs.io/en/latest/index.html)
* Doublet detection tools:
  * [Scds script](demuxafy_scds_runner.qsub) using [scds](https://github.com/kostkalab/scds)
  * [ScDblFinder script](demuxafy_scdblfinder_runner.qsub) uses [scDblFinder](https://github.com/plger/scDblFinder)
* [Combiner script](demuxafy_combiner.qsub) combines all results (atm: vireo, scds, scDblFinder) using majority voting

### Pre-processing

* [preprocessing scripts folder](preprocessing/)

### Post-processing

Update maxi pool vireo results as correct ones, and non-maxi pool ones as outdated.
This only affects S0041, S0043, S0044, S0045, S0046, S0047.

```bash
cd /directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/vireo_output_no_cb/

# move single pools to old
mv S0041/ old_S0041
mv S0043/ old_S0043
mv S0044/ old_S0044
mv S0045/ old_S0045
mv S0046/ old_S0046

# move maxi pools to current
mv maxi_pool_S0041/ S0041/
mv maxi_pool_S0043/ S0043/
mv maxi_pool_S0044/ S0044/
mv maxi_pool_S0045/ S0045/
mv maxi_pool_S0046/ S0046/
mv maxi_pool_S0047/ S0047/
```


<!---
* [demuxafy_dropulation_manual_summary.R](demuxafy_dropulation_manual_summary.R)

### other tools

* [Demuxalot script](demuxafy_demuxalot_runner.qsub) using [Demuxalot](https://pypi.org/project/demuxalot/), requires genotypes
* [Dropulation script](demuxafy_dropulation_runner.qsub) using [Dropulation](https://github.com/broadinstitute/Drop-seq/blob/master/doc/Census-seq_Computational_Protcools.pdf), requires genotypes
* [DoubletDetection script](demuxafy_doubletdetection_runner.qsub) uses [DoubletDetection](https://github.com/JonathanShor/DoubletDetection)
--->
