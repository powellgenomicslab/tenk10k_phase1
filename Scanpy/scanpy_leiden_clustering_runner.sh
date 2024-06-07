#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=800G
#$ -l tmp_requested=40G
#$ -N scanpy_leiden
#$ -q long.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/logs/scanpy_leiden.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/logs/scanpy_leiden.stdout
#$ -m ae
#$ -M b.bowen@garvan.org.au

PY_PATH="/directflow/SCCGGroupShare/projects/blabow/.conda/envs/cellregmap_notebook_bb/bin"
PYSCRIPT="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/Scanpy/scanpy_leiden_clustering.py"

cd ${TMPDIR}
# Run python script
echo "${PY_PATH}/python ${PYSCRIPT}"

${PY_PATH}/python ${PYSCRIPT}