#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=20G
#$ -l tmp_requested=40G
#$ -pe smp 12
#$ -N coloc_ukbb
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc_ukbb.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc_ukbb.stdout
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

# RUN COLOC FOR UKBB BLOOD TRAITS ----

GWAS_DIR_UKBB=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/gwas/gymrek-ukbb-snp-str-gwas-catalogs/chr-specific/

OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/coloc_results/${CELLTYPE}/blood_traits/
mkdir -p ${OUTDIR}

for GWAS_FILE_PATH in ${GWAS_DIR_UKBB}/*; do

    gwas=$(basename "${GWAS_FILE_PATH}")
    GWAS_PHENO_NAME=$(echo $gwas | sed 's/_snp_str_gwas_results_hg38_chr[0-9]\+\.tab\.gz$//') # parse phenotype name from the file name 

    parallel -j 12 --verbose run_coloc ${GWAS_FILE_PATH} ${GWAS_PHENO_NAME} {} ${CELLTYPE} ${OUTDIR} ::: $(ls /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/plink_fixed_pheno/${CELLTYPE}/*.glm.linear)

done
