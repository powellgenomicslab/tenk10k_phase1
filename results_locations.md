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

### Demuxafy

* scds results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/scds_output/{SAMPLE}/```
* scDblFinder results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/scdblfinder_output/{SAMPLE}/```
* vireo results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/vireo_output_no_cb/{SAMPLE}/```
* demuxafy combined results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/demuxafy/combined_output_scds_scdblfinder_vireo_no_cb/{SAMPLE}/```

### Other

* Consortium WG2 Azimuth results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/{SAMPLE}/step2_azimuth/```
* Consortium WG2 hierarchical scPred results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scpred/{SAMPLE}/step3_hierscpred/```
* Celltypist results: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/celltypist/```
* Scanpy / AnnData objects with metadata: ```/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/scanpy_objects_w_metadata/```
