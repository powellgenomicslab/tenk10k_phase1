#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=20G
#$ -l tmp_requested=40G
#$ -pe smp 12
#$ -N coloc_disease
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc_disease.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc_disease.stdout
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

RUN COLOC FOR DISEASE TRAITS ----

GWAS_DIR_DISEASE=/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/gwas/

pheno_names=(
    "alzheimer_GCST90027158"
    "breastca_GCST004988"
    "colorectalca_GCST90129505"
    "covid_GCST011071"
    "ibd_liu2023"
    "NHL_GCST90011819"
    "lungca_GCST004748"
    "lymphoma_GCST90018878"
    "parkinson_GCST009325"
    "prostateca_GCST90274713"
    "ra_GCST90132223"
    "sle_GCST003156"
    "myeloproliferative_GCST90000032"
    "lymphocytic_leukemia_GCST90011814"
    "nephrotic_GCST90258619"
    "kiryluk_IgAN"
)

pheno_files=(
    "${GWAS_DIR_DISEASE}/GCST90027158.h_parsed.tsv"
    "${GWAS_DIR_DISEASE}/GCST004988.h_parsed.tsv"
    "${GWAS_DIR_DISEASE}/colorectalca_GCST90129505_parsed.tsv"
    "${GWAS_DIR_DISEASE}/GCST011071_parsed.tsv"
    "${GWAS_DIR_DISEASE}/ibd_EAS_EUR_SiKJEF_meta_IBD.tsv"
    "${GWAS_DIR_DISEASE}/NHL_GCST90011819_parsed.tsv"
    "${GWAS_DIR_DISEASE}/GCST004748.h_parsed.tsv"
    "${GWAS_DIR_DISEASE}/GCST90018878.h_parsed.tsv"
    "${GWAS_DIR_DISEASE}/GCST009325.h_parsed.tsv"
    "${GWAS_DIR_DISEASE}/GCST90274713.h_parsed.tsv"
    "${GWAS_DIR_DISEASE}/GCST90132223_parsed.tsv"
    "${GWAS_DIR_DISEASE}/bentham_2015_26502338_sle_parsed.tsv"
    "${GWAS_DIR_DISEASE}/myeloproliferative_GCST90000032_parsed.tsv"
    "${GWAS_DIR_DISEASE}/lymphocytic_leukemia_GCST90011814_parsed.tsv"
    "${GWAS_DIR_DISEASE}/nephrotic_GCST90258619_parsed.tsv"
    "${GWAS_DIR_DISEASE}/Kiryluk_IgAN_parsed.tsv"
)

OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/${CELLTYPE}/disease_traits/
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