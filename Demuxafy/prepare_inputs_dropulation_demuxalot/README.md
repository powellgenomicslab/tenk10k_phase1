# Additional Data Preparation

Demuxalot and Dropulation (the two demultiplexing tools used here) need genotype files (as well as the CellRanger outputs), which need a little more preparation.

## Expression files

```BAM``` files and ```barcodes.tsv``` are outputted directly via CellRanger.

## Genotypes

### Individual files

First, get the list of individuals from the ```VCF``` using ```bcftools```:

```bash
module use /share/ClusterShare/apps/brenner/Modules/modulefiles
module load bcftools
bcftools query -l /directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/genotype_hg38/Merged_MAF0.05_hg38_nochr.vcf.gz > /share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/donor_list.txt
```

### Restricting VCF

I may potentially need to reduce the ```VCF``` file to only exonic variants following the documentation: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DataPrep.html#filter-for-snps-overlapping-exons

## Gene annotation (Dropulation only)

Downloaded the "Basic gene annotation" ```GTF``` from https://www.gencodegenes.org/human/ -- it is definitely the correct version (v44), potentially not quite exactly the same file used in the alignment (there are many files in there, plus I think CellRanger may modify the GTF??).

