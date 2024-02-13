# CellBender for ambient RNA detection

We run CellBender ([paper](https://www.nature.com/articles/s41592-023-01943-7), [docs](https://cellbender.readthedocs.io/en/latest/index.html)) to detect the levels of ambient RNA present in our droplets, and to estimate empty vs non-empty droplets.

We consider default parameters except for a decreased ```learning rate``` value (1e-6 instead of the default 1e-4).

At present we do not perform any hard filters based on the Cellbender results, but annotate our cells based on the % ambient RNA identified for future reference.
We do this as part of the [script]() where we annotate each library's AnnData object with results from different tools.

Note that we are running v0.3.0 upon recommendation of the lead author, Stephen Fleming (there were some issues with v0.3.1).

## Brief description of scripts included here

* [cellbender_runner.qsub](cellbender_runner.qsub): qsub script to submit job for running cellbender across all samples
