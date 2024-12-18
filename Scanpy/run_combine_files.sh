#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=1500G
#$ -l tmp_requested=40G
#$ -N combine_scanpy
#$ -q long.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/logs/combine_files.err
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/logs/combine_files.log
#$ -m ae
#$ -M b.bowen@garvan.org.au

PY_PATH="/directflow/SCCGGroupShare/projects/blabow/.conda/envs/cellregmap_notebook_bb/bin"

PYSCRIPT="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/Scanpy/combine_files_add_gene_info.py"

## SCRIPT TO FILTER COMBINED OBJECT
# PYSCRIPT="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/Scanpy/filter_combined_scanpy.py"

echo "${PY_PATH}/python ${PYSCRIPT}"

${PY_PATH}/python ${PYSCRIPT}