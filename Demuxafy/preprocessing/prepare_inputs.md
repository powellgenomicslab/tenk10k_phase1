# Simple pre-processing steps

## Extract cellranger libraries 

This is used for scds, scDblFinder to know what samples to conisder when parallelising.

```bash
cd /share/ScratchGeneral/anncuo/tenk10k/data_processing/
ls /directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231013_tenk10k_gencode44/cellranger_outs/ > cellranger_outs_231013.txt
ls /directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231213_tenk10k_gencode44/cellranger_outs/ > cellranger_outs_231213.txt
ls /directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231214_tenk10k_gencode44/cellranger_outs/ > cellranger_outs_231214.txt
```

## VCF manipulation

### add "chr" to VCF to match BAM files

First, copy to own folder and unzip:

```bash
cp /directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/genotype_hg38/Merged_MAF0.05_hg38_nochr.vcf.gz /share/ScratchGeneral/anncuo/tenk10k/data_processing/genotypes/Merged_MAF0.05_hg38_nochr.vcf.gz
gunzip /share/ScratchGeneral/anncuo/tenk10k/data_processing/genotypes/Merged_MAF0.05_hg38_nochr.vcf.gz /share/ScratchGeneral/anncuo/tenk10k/data_processing/genotypes/Merged_MAF0.05_hg38_nochr.vcf
```

Then, add "chr":

```bash
awk '{
        if($0 !~ /^#/)
            print "chr"$0;
        else if(match($0,/(##contig=<ID=)(.*)/,m))
            print m[1]"chr"m[2];
        else print $0
      }' Merged_MAF0.05_hg38_nochr.vcf > Merged_MAF0.05_hg38_chr.vcf
```

### zip and index

```bash
module use /share/ClusterShare/apps/brenner/Modules/modulefiles
module load bcftools
module load htslib
bgzip -c Merged_MAF0.05_hg38_chr.vcf > Merged_MAF0.05_hg38_chr.vcf.gz
tabix -p vcf Merged_MAF0.05_hg38_chr.vcf.gz
```


### Individual files directly from VCF

First, get the list of individuals from the ```VCF``` using ```bcftools```:

```bash
module use /share/ClusterShare/apps/brenner/Modules/modulefiles
module load bcftools
bcftools query -l /directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/genotype_hg38/Merged_MAF0.05_hg38_nochr.vcf.gz > /share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/donor_list.txt
```

<!--- 
### Restricting VCF

I may potentially need to reduce the ```VCF``` file to only exonic variants following the documentation: https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/DataPrep.html#filter-for-snps-overlapping-exons

## Gene annotation (Dropulation only)

Downloaded the "Basic gene annotation" ```GTF``` from https://www.gencodegenes.org/human/ -- it is definitely the correct version (v44), potentially not quite exactly the same file used in the alignment (there are many files in there, plus I think CellRanger may modify the GTF??). --->

