#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=400G
#$ -l tmp_requested=40G
#$ -N expr_pcs
#$ -q long.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA_expr_pcs.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA_GeNA_expr_pcs.stdout
#$ -t 4-4
#$ -m ae
#$ -M b.bowen@garvan.org.au

# NOTE: some of the cell types require more memory than others

PY_PATH="/directflow/SCCGGroupShare/projects/blabow/.conda/envs/gena-env/bin"
PYSCRIPT=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/04.5_make_expression_pcs.py

i=${SGE_TASK_ID};
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/unique_cell_types_wg2_scpred.txt)
LOG="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/${CELLTYPE}_${SGE_TASK_ID}_expr_pcs.log"

# Do the main job
${PY_PATH}/python ${PYSCRIPT} ${CELLTYPE} &>> ${LOG}
