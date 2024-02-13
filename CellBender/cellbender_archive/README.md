## Brief description of scripts included here

* [cellbender_runner.qsub](cellbender_runner.qsub): qsub script to submit job for running cellbender across all samples
* comparison scripts to evaluate cellbender results
  * [celltypist_confidence_cells_uncorrected.py](celltypist_confidence_cells_uncorrected.py): comparison of average cell typing confidence score using CellTypist, when comparing cells estimated by cellranger vs cellbender (**uncorrected counts**), so evaluating solely differences in cell inclusion criteria
  * [celltypist_confidence_cells_corrected.py](celltypist_confidence_cells_corrected.py): same but for **corrected counts**, evaluating both inclusion criteria of cells, and effects of ambient RNA correction
  * [cellbender_celltypist_confidence_plot_across_samples.R](cellbender_celltypist_confidence_plot_across_samples.R): plot results from above scripts
* [vireo_runner_cellbender_file.qsub](vireo_runner_cellbender_file.qsub): qsub script to submit job for running vireo considering all cells identified as such by cellbender (not yet tested)
* [cellbender_celltypist_confidence_plot_across_samples.R](cellbender_celltypist_confidence_plot_across_samples.R): compare vireo results when using cellranger vs cellbender barcodes
