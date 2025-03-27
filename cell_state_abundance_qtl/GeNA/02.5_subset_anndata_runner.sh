#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=1000G
#$ -l tmp_requested=40G
#$ -N subset_anndata
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA_prep_data.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA_prep_data.stdout
#$ -t 1-8
#$ -m ae
#$ -M b.bowen@garvan.org.au

# NOTE: some of the cell types require more memory than others

PY_PATH="/directflow/SCCGGroupShare/projects/blabow/.conda/envs/gena-env/bin"
PYSCRIPT="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/02.5_subset_anndata.py"


i=${SGE_TASK_ID};

RESOLUTION="major_cell_types"
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)

# RESOLUTION="minor_cell_types"
# CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/unique_cell_types_wg2_scpred.txt)

LOG="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/${RESOLUTION}_${CELLTYPE}_${SGE_TASK_ID}_subset_anndata.log"

${PY_PATH}/python ${PYSCRIPT} ${CELLTYPE} ${RESOLUTION} &>> ${LOG}

#2,3,4,8

# Note:
# 
# CD4 T cells and CD8 T cells require very high memory to run this script;
# I think this is due to the all celltype anndata persisting in memory after being subsetted 
# can likely reduce this and 
# 