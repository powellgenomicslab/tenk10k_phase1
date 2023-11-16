# Additional Data Preparation

Demuxalot and Dropulation (the two demultiplexing tools used here) need genotype files (as well as the CellRanger outputs), which need a little more preparation.

BAM files and barcodes are outputted directly via CellRanger.

## Individual files

First, get the list of individuals from the VCF using bcftools:

```bash
module use /share/ClusterShare/apps/brenner/Modules/modulefiles
module load bcftools
bcftools query -l /directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/genotype_hg38/Merged_MAF0.05_hg38_nochr.vcf.gz > /share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/donor_list.txt
```

## Restricting VCF

I may potentially need to reduce the VCF file to only exonic variants following the documentation: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DataPrep.html#filter-for-snps-overlapping-exons
