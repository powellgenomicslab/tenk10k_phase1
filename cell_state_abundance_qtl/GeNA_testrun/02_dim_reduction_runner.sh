#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=800G
#$ -l tmp_requested=40G
#$ -N mad_dim_reduc
#$ -q long.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/mad_dim_reduc.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/mad_dim_reduc.stdout
#$ -m ae
#$ -M b.bowen@garvan.org.au

PY_PATH="/directflow/SCCGGroupShare/projects/blabow/.conda/envs/gena-env/bin"
PYSCRIPT="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/02_dim_reduction.py"

# Run python script
echo "${PY_PATH}/python3 ${PYSCRIPT}"

${PY_PATH}/python ${PYSCRIPT}