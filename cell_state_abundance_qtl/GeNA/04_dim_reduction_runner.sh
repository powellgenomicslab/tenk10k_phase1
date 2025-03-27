#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=800G
#$ -l tmp_requested=40G
#$ -N GeNA_dimreduc
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA_dim_reduction.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA_dim_reduction.stdout
#$ -t 9-9
#$ -m ae
#$ -M b.bowen@garvan.org.au

# NOTE: some of the cell types require more memory than others

PY_PATH="/directflow/SCCGGroupShare/projects/blabow/.conda/envs/gena-env/bin"
PYSCRIPT="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/04_dim_reduction.py"

i=${SGE_TASK_ID};
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)
RESOLUTION="major_cell_types"

LOG="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/${RESOLUTION}_${CELLTYPE}_${SGE_TASK_ID}_dim_reduction.log"

# Do the main job
${PY_PATH}/python ${PYSCRIPT} ${CELLTYPE} ${RESOLUTION} &>> ${LOG}
