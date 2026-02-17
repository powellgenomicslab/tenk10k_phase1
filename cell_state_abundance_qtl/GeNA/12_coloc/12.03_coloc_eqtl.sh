#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=20G
#$ -l tmp_requested=40G
#$ -pe smp 12
#$ -N coloc
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc_eqtl.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc_eqtl.stdout
#$ -t 1-7
#$ -m ae
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate coloc-env

CELLTYPE=$(sed "${SGE_TASK_ID}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)

# source coloc runner bash function 
source /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/12_coloc/12.04_coloc_function.sh
export -f run_coloc_eqtl

# RUN COLOC for eQTL GWAS ----
# coloc is run for all csaQTL for a major cell type vs eQTL for minor cell types within that major cell type

# EQTL_DIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/from_smr/reformatted_for_coloc
EGENES_FILE_PATH=/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/saige_qtl/december24_freeze/
CT_TSV=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/plotting_notebooks/overview_figures/manuscript_figures/colour_palette_table.tsv
# use awk to get minor celltypes for each major cell type. filter the above table so column 10 matches $CELLTYPE, return column 1 as an array
# MINOR_CELLTYPES=$(awk -v ct="$CELLTYPE" 'BEGIN{FS="\t"} $10==ct {print $1}' $CT_TSV)

# Its quite fast so just run for all minor cell types
MINOR_CELLTYPES=$(awk 'BEGIN{FS="\t"} {print $1}' $CT_TSV)

OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/${CELLTYPE}/eqtl/
mkdir -p ${OUTDIR}

cmd="parallel -j 12 --verbose run_coloc_eqtl ${EGENES_FILE_PATH} {1} ${CELLTYPE} {2} ${OUTDIR} ::: $MINOR_CELLTYPES ::: $(ls /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/plink_fixed_pheno/${CELLTYPE}/*.glm.linear)"
echo $cmd 
eval $cmd
