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

# RUN COLOC for Vukovic et al UKBB GWAS ----

GWAS_DIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/from_smr/reformatted_for_coloc

pheno_names=(
    Basophil_count
    eosinophil_count
    hematocrit
    hemoglobin_concentration
    hlr_reticulocyte_count
    hlr_reticulocyte_percentage_of_red_cells
    immature_fraction_of_reticulocytes
    lymphocyte_count
    mean_corpuscular_hemoglobin_concentration
    mean_corpuscular_hemoglobin
    mean_corpuscular_volume
    monocyte_count
    monocyte_percentage_of_white_cells
    mean_platelet_volume
    mean_reticulocyte_volume
    mean_sphered_corpuscular_volume
    neutrophil_count
    neutrophil_percentage_of_white_cells
    plateletcrit
    platelet_distribution_width
    platelet_count
    red_blood_cell_count
    red_cell_distribution_width
    reticulocyte_count
    reticulocyte_fraction_of_red_cells
    white_blood_cell_count
)

pheno_files=(
    "${GWAS_DIR}/baso.ma"
    "${GWAS_DIR}/eo.ma"
    "${GWAS_DIR}/hct.ma"
    "${GWAS_DIR}/hgb.ma"
    "${GWAS_DIR}/hlr.ma"
    "${GWAS_DIR}/hlr_p.ma"
    "${GWAS_DIR}/irf.ma"
    "${GWAS_DIR}/lymph.ma"
    "${GWAS_DIR}/mchc.ma"
    "${GWAS_DIR}/mch.ma"
    "${GWAS_DIR}/mcv.ma"
    "${GWAS_DIR}/mono.ma"
    "${GWAS_DIR}/mono_p.ma"
    "${GWAS_DIR}/mpv.ma"
    "${GWAS_DIR}/mrv.ma"
    "${GWAS_DIR}/mscv.ma"
    "${GWAS_DIR}/neut.ma"
    "${GWAS_DIR}/neut_p.ma"
    "${GWAS_DIR}/pct.ma"
    "${GWAS_DIR}/pdw.ma"
    "${GWAS_DIR}/plt.ma"
    "${GWAS_DIR}/rbc.ma"
    "${GWAS_DIR}/rdw_cv.ma"
    "${GWAS_DIR}/ret.ma"
    "${GWAS_DIR}/ret_p.ma"
    "${GWAS_DIR}/wbc.ma"
)

OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/${CELLTYPE}/extra_traits/
mkdir -p ${OUTDIR}

i=0
for GWAS_FILE_PATH in "${pheno_files[@]}"; do
    gwas=$(basename "${GWAS_FILE_PATH}")
    GWAS_PHENO_NAME="${pheno_names[$i]}"
    cmd="parallel -j 12 --verbose run_coloc ${GWAS_FILE_PATH} ${GWAS_PHENO_NAME} {} ${CELLTYPE} ${OUTDIR} ::: $(ls /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/plink_fixed_pheno/${CELLTYPE}/*.glm.linear)"
    eval $cmd
    ((i+=1))
done