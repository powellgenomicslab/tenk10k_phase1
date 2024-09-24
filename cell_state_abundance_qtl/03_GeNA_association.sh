#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=250G
#$ -l tmp_requested=40G
#$ -N GeNA
#$ -q long.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA.stdout
#$ -m ae
#$ -M b.bowen@garvan.org.au

# conda environment 
. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate gena-env

# go to GeNA repo
cd /directflow/SCCGGroupShare/projects/blabow/software/GeNA

./GeNA.sh -s /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/scDataObject.dimreduc.h5ad \
    -b 'False' \
    -g /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/Merged_MAF0.05_hg38_chr_snps \
    -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/ \
    -c 'age,sex,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5'

# TODO:
# [] add sex back in once I replace missing sex as zero 
# [] get sample-level PC's (either take avg of harmony PC's per-sample OR do PCA on the pseudobulk expression matrix...) 

# NOTE: removed sex as covariate for testing if the NA value is what's causing the error, it ran OK when I removed sex
            # therefore probably was erroring because of missing sex

#    -c 'age,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6,geno_PC7,geno_PC8,geno_PC9,geno_PC10,geno_PC11,geno_PC12,geno_PC13,geno_PC14,geno_PC15,geno_PC16'

# TODO:
#   Currently not possible to use batch correction in GeNA
#   [x] figure out how to add the batch correction, given the many-to-many relationship between batches and samples (the same sample is present in multiple seq libraries...) 
#       * look at GeNA-applied github for potential answers, or ask author

# USAGE:
#   -s the path to the single-cell object in MultiAnnData format as demonstrated above
#   -b whether or not to use sample assignments to batches stored in a d.samplem column titled 'batch' to remove batchy neighborhoods
#   -g a path to genotyping data in PLINK2 format
#   -o a path to an output folder for GeNA sumstats
#   -c a list of covariates to control for when running GeNA. Expected format is a comma-delimited list of column names stored in d.samplem
#   -k (Optional) a path to a file with user-specified values of k for GeNA to consider, formatted with one value of k on each line
