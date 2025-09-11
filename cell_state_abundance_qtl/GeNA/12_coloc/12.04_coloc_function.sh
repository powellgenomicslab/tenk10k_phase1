#! /bin/bash

# function to run analysis for 1 csaQTL SNP / GWAS 
# to be run with parallel 

run_coloc() {
    local gwas_file_path=$1
    local gwas_pheo_name=$2 
    local csaqtl_file_path=$3
    local CELLTYPE=$4
    local OUTDIR=$5

    python /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/12_coloc/12.05_coloc.py \
        --gwas-file-path ${gwas_file_path} \
        --gwas-pheno-name ${gwas_pheo_name} \
        --csaqtl-file-path ${csaqtl_file_path} \
        --celltype ${CELLTYPE} \
        --outdir ${OUTDIR} 
}

# function to run analysis for 1 csaQTL SNP / eQTL for one major cell type / minor cell type combination 
run_coloc_eqtl() {
    local egenes_file_path=$1
    local celltype_eqtl=$2 
    local celltype_csaqtl=$3
    local csaqtl_file_path=$4
    local OUTDIR=$5

    cmd="python /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/12_coloc/12.05_coloc_eqtl.py \
            --egenes-files-path ${egenes_file_path} \
            --celltype-eqtl ${celltype_eqtl} \
            --celltype-csaqtl ${celltype_csaqtl} \
            --csaqtl-file-path ${csaqtl_file_path} \
            --outdir ${OUTDIR}"

    echo $cmd
    eval $cmd
}