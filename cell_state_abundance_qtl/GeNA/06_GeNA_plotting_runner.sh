#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=120G
#$ -N plotting
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA_plotting.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/GeNA_plotting.stdout
#$ -t 1-8
#$ -m ae
#$ -M b.bowen@garvan.org.au

# -----------
# Parameters:
# -----------

# ANALYSIS_NAME="no_expr_pc_covars"
ANALYSIS_NAME="with_expr_pc_covars"
RESOLUTION="major_cell_types"

# -----------

i=${SGE_TASK_ID};
# CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/unique_cell_types_wg2_scpred.txt)
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate r-python-dev

Rscript /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/06_GeNA_plotting.R ${RESOLUTION} ${CELLTYPE} ${ANALYSIS_NAME}