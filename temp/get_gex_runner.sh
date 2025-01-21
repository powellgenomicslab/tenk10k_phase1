#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=1500G
#$ -l tmp_requested=40G
#$ -N get_gex
#$ -q long.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/get_gex.err
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/get_gex.log
#$ -m ae
#$ -M b.bowen@garvan.org.au

PY_PATH="/directflow/SCCGGroupShare/projects/blabow/.conda/envs/cellregmap_notebook_bb/bin"

# first pass harmony + UMAP
PYSCRIPT="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/get_gene_expression.py"

cd ${TMPDIR}
# Run python script
echo "${PY_PATH}/python ${PYSCRIPT}"

${PY_PATH}/python ${PYSCRIPT}