#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=350G
#$ -N saveumap
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/saveumap.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/saveumap.stdout
#$ -t 1-8
#$ -m ae
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate gena-env

# -----------
# Parameters:
# -----------

# these args are just used to put re-runs of the analysis in different output directories 
RESOLUTION="major_cell_types"

i=${SGE_TASK_ID};
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)

# run python script to calculate sample and neighbourhood level phenotypes 
python /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/04.6_save_umaps.py ${CELLTYPE} ${RESOLUTION} ${ANALYSIS_NAME} ${COVS}

