#! /bin/bash

# function to run analysis for 1 SNP / GWAS 
run_coloc() {
    local gwas_file_path=$1
    local gwas_pheo_name=$2 
    local csaqtl_file_path=$3
    local CELLTYPE=$4
    local OUTDIR=$5

    python /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/12_coloc/12.2_coloc.py \
        --gwas-file-path ${gwas_file_path} \
        --gwas-pheno-name ${gwas_pheo_name} \
        --csaqtl-file-path ${csaqtl_file_path} \
        --celltype ${CELLTYPE} \
        --outdir ${OUTDIR} 
}