#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=350G
#$ -N get_pheno
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/get_pheno.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/get_pheno.stdout
#$ -t 1-8
#$ -m ae
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate gena-env

# -----------
# Parameters:
# -----------

# these args are just used to put re-runs of the analysis in different output directories 
ANALYSIS_NAME="no_expr_pc_covars"
RESOLUTION="major_cell_types"
# covariates used in the gena model 
COVS="sex,age,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6,geno_PC7,BioHEART"

i=${SGE_TASK_ID};
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)

# run python script to calculate sample and neighbourhood level phenotypes 
python /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/08.1_get_pheno.py ${CELLTYPE} ${RESOLUTION} ${ANALYSIS_NAME} ${COVS}

