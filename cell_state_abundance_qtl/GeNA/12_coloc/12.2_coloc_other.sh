#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=20G
#$ -l tmp_requested=40G
#$ -pe smp 12
#$ -N coloc
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc.stdout
#$ -t 1-7
#$ -m ae
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate coloc-env

ANALYSIS_NAME="no_expr_pc_covars"
RESOLUTION="major_cell_types"

CELLTYPE=$(sed "${SGE_TASK_ID}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)

# source coloc runner bash function 
source /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/12_coloc/coloc_function.sh
export -f run_coloc

# RUN COLOC FOR DISEASE TRAITS ----

GWAS_DIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/from_smr/reformatted_for_coloc

pheno_names=(
    "mono_p"
    "mono"
)

pheno_files=(
    "${GWAS_DIR}/mon_P.ma"
    "${GWAS_DIR}/mono.ma"
)

OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/${CELLTYPE}/extra_traits/
mkdir -p ${OUTDIR}

i=0
for GWAS_FILE_PATH in "${pheno_files[@]}"; do
    gwas=$(basename "${GWAS_FILE_PATH}")
    GWAS_PHENO_NAME="${pheno_names[$i]}"
    # echo $i
    # echo $GWAS_PHENO_NAME
    # echo $GWAS_FILE_PATH
    cmd="parallel -j 12 --verbose run_coloc ${GWAS_FILE_PATH} ${GWAS_PHENO_NAME} {} ${CELLTYPE} ${OUTDIR} ::: $(ls /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/plink_fixed_pheno/${CELLTYPE}/*.glm.linear)"
    # echo $cmd
    eval $cmd
    ((i+=1))
done