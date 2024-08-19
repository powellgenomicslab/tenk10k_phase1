# Data processing results location

Keeping track of where the results are currently (on Brenner).

The goal is to eventually move everything to a new directory at ```/directflow/SCCGGroupShare/projects/tenk10k_phase1```.

## CellBender results

  * Anna: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/cellbender_output_smaller_learning_rate/{SAMPLE}/```
  * Blake: ```/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cellbender/output/smaller_learning_rate/{SAMPLE}/```

Files in dir:

* ```cellbender_output_cell_barcodes.csv```
* ```cellbender_output.log```
* ```cellbender_output_posterior.h5```
* ```cellbender_output_filtered.h5```
* ```cellbender_output_metrics.csv```
* ```cellbender_output_report.html```
* ```cellbender_output.h5```
* ```cellbender_output.pdf```

Actually used:

* ```cellbender_output_cell_barcodes.csv``` at least used to be used for vireo?
* ```cellbender_output_metrics.csv``` for plotting
* one of the ```.h5``` files?

TO DO: clarify and zip the rest

### Scripts pointing to CellBender results

Scripts that point to those locations (and therefore need to be changed):

* [CellBender runner](CellBender/cellbender_runner.qsub)
* [Scanpy adding info to AnnData](Scanpy/add_metadata_per_sample_no_norm.py)

## Demuxafy
### scds results:

* Anna: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/scds_output/{SAMPLE}/```
* Blake: ```/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/scds/output/${SAMPLE}```
Files in dir: 
* `scds_doublets_singlets.tsv`
* `scds_doublet_summary.tsv`
### scDblFinder results:
*  Anna: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/scdblfinder_output/{SAMPLE}/```
* Blake: ```/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/scdblfinder/output/${SAMPLE}```
Files in dir:
* `scDblFinder_doublets_singlets.tsv` 
* `scDblFinder_doublet_summary.tsv`

### vireo results:
* Anna:  ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/vireo_output_no_cb/{SAMPLE}/```
* Blake: `/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/vireo_all_samples/output/${SAMPLE}`
Files in dir:
* `cellSNP.base.vcf.gz`
* `cellSNP.tag.DP.mtx`
* `donor_subset.vcf`
* `prob_doublet.tsv.gz`
* `summary.tsv`
* `cellSNP.samples.tsv` 
* `cellSNP.tag.OTH.mtx` 
* `fig_GT_distance_estimated.pdf`
* `prob_singlet.tsv.gz`
* `cellSNP.tag.AD.mtx` 
* `donor_ids.tsv`       
* `_log.txt`
* `prop_ambient.tsv`
### demuxafy combined results: 
* Anna:  ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/combined_output_scds_scdblfinder_vireo_no_cb/{SAMPLE}/```
* Blake: `/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/combine_results/output/combined_output_scds_scdblfinder_vireo_no_cb/${SAMPLE}`
Files in dir: 
* `combined_results_assignment_summary.tsv`
* `combined_results_Singlets_upset.pdf`
* `combined_results_demultiplexing_summary.tsv`
* `combined_results_summary.tsv`
* `combined_results_droplet_type_summary.tsv `
* `combined_results.tsv`
* `combined_results_Singlets_upset_donor_assignment.pdf `
* `combined_results_w_combined_assignments.tsv`
* `combined_results_Singlets_upset_droplet_type.pdf`
## Celltyping

### Unfiltered Seurat Objects: 
* Anna:
* Blake: `/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scpred/seurat_objects_unfiltered/`

### scPred

#### Consortium WG2 Azimuth results:
 * Anna:  ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/{SAMPLE}/step2_azimuth/```
  * Blake: `/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scpred/${SAMPLE}/step2_azimuth/`
  Files in dir:
  * `azimuth.RDS`
  * `azimuth_ref_spca.png`
  * `azimuth_ref_umap.png`
#### Consortium WG2 hierarchical scPred results:
 * Anna:  ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/{SAMPLE}/step3_hierscpred/```
 * Blake: `/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scpred/${SAMPLE}/step3_hierscpred/`
Files in dir:
* hier_scpred.RDS
### Combined scpred + azimuth results

* Anna: 
* Blake: `/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scpred/${SAMPLE}"`
Files in dir: 
* `combined_metadata.csv`
### Celltypist

#### Celltypist results:
 * Anna: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/celltypist/```
 * Blake: `/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/celltypist/output/`
files in dir:
* `${SAMPLE}_celltypist_predicted.h5`
### Scanpy

#### Scanpy / AnnData objects with metadata:
 * Anna:  ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scanpy_objects_w_metadata/```
 * Blake: `/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/scanpy_objects_w_metadata/${SAMPLE}`
Files in dir: 
* `${SAMPLE}_w_metadata_donor_info.h5ad`
