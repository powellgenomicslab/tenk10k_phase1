#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=50G
#$ -l tmp_requested=4G
#$ -N perm_GeNA
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA.stdout
#$ -m ae
#$ -M b.bowen@garvan.org.au
#$ -t 1-1000
#$ -tc 250

# job array across cell types 
CELLTYPE="Monocyte"
RESOLUTION="major_cell_types"

# parameters 
OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}
SCDATA=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/${RESOLUTION}/${CELLTYPE}_scDataObject.dimreduc.pca.h5ad # THESE HAVE PSEUDOBULK PCA DIMs - currently only works with these because the other version does not have Cohort in the metadata
# SCDATA=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/${RESOLUTION}/${CELLTYPE}_scDataObject.dimreduc.h5ad # WITHOUT PSEUDOBULK PCs, which I decided not to use anyway, note to remove pseudobulking step, add cohort mappings to single cell objects in earlier script 

GENOTYPES=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15

# conda environment 
. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate gena-env

# go to GeNA repo
cd /directflow/SCCGGroupShare/projects/blabow/software/GeNA

###############################################################
# Run a bunch of permutations for a single variant of interest 
###############################################################

VARIANT="15:39687137:C:T"

mkdir -p ${OUTDIR}/${CELLTYPE}/single_variant_permutations/

i=${SGE_TASK_ID};

PERMUTED_GTYPE_SINGLE_VARIANT=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted/variant_permutations/${VARIANT}_perm${i}

echo " Runing GeNA for ${VARIANT} genotype permutation # ${i}"

# getting an error that output directory doesn't exist
mkdir -p ${OUTDIR}/${CELLTYPE}/single_variant_permutations/single_variant_${VARIANT}_perm${i}

./GeNA.sh -s ${SCDATA} \
    -b 'False' \
    -g ${PERMUTED_GTYPE_SINGLE_VARIANT} \
    -o ${OUTDIR}/${CELLTYPE}/single_variant_permutations/single_variant_${VARIANT}_perm${i}/ \
    -k /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted/variant_permutations/k.txt \
    -c 'sex,age,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6,geno_PC7,BioHEART' # covariates 

