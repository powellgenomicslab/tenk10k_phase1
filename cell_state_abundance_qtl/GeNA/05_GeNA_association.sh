#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=400G
#$ -l tmp_requested=40G
#$ -N GeNA
#$ -q long.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA.stdout
#$ -m ae
#$ -M b.bowen@garvan.org.au
#$ -t 7-7

# job array across cell types 
i=${SGE_TASK_ID};
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)
RESOLUTION="major_cell_types"

LOG=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/${RESOLUTION}_${CELLTYPE}_${SGE_TASK_ID}_GeNA_association.log

# parameters 
OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}
SCDATA=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/${RESOLUTION}/${CELLTYPE}_scDataObject.dimreduc.pca.h5ad # THESE HAVE PSEUDOBULK PCA DIMs - currently only works with these because the other version does not have Cohort in the metadata
# SCDATA=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/${RESOLUTION}/${CELLTYPE}_scDataObject.dimreduc.h5ad # WITHOUT PSEUDOBULK PCs, which I decided not to use anyway, note to remove pseudobulking step, add cohort mappings to single cell objects in earlier script 

GENOTYPES=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15
PERMUTED_GTYPES=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted/merged_common_variants_standard_chr_geno_0.15_singlecell_cohort_perm1
# conda environment 
. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate gena-env

# make outdirs
mkdir -p ${OUTDIR}/${CELLTYPE}/no_expr_pc_covars/ 
mkdir -p ${OUTDIR}/${CELLTYPE}/with_expr_pc_covars/ 
mkdir -p ${OUTDIR}/${CELLTYPE}/no_expr_pc_covars_perm/
mkdir -p ${OUTDIR}/${CELLTYPE}/with_expr_pc_covars_perm/

# go to GeNA repo
cd /directflow/SCCGGroupShare/projects/blabow/software/GeNA

(
    echo "Running GeNA for ${CELLTYPE} - NO expression PC covars"
    echo ""

    GeNA association: NO expression PC covariates

    run association 
    ./GeNA.sh -s ${SCDATA} \
        -b 'False' \
        -g ${GENOTYPES} \
        -o ${OUTDIR}/${CELLTYPE}/no_expr_pc_covars/ \
        -c 'sex,age,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6,geno_PC7,BioHEART' # covariates 

    # run association with random genotypes
    ./GeNA.sh -s ${SCDATA} \
        -b 'False' \
        -g ${PERMUTED_GTYPES} \
        -o ${OUTDIR}/${CELLTYPE}/no_expr_pc_covars_perm/ \
        -c 'sex,age,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6,geno_PC7,BioHEART' # covariates 



) &> ${LOG}


# USAGE:
#   -s the path to the single-cell object in MultiAnnData format as demonstrated above
#   -b whether or not to use sample assignments to batches stored in a d.samplem column titled 'batch' to remove batchy neighborhoods
#   -g a path to genotyping data in PLINK2 format
#   -o a path to an output folder for GeNA sumstats
#   -c a list of covariates to control for when running GeNA. Expected format is a comma-delimited list of column names stored in d.samplem
#   -k (Optional) a path to a file with user-specified values of k for GeNA to consider, formatted with one value of k on each line


