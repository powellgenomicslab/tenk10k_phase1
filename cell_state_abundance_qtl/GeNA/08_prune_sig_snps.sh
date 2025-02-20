#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=10G
#$ -N prune
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/prune_sig_snps.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/prune_sig_snps.stdout
#$ -t 1-8
#$ -pe smp 10
#$ -m ae
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate gena-env

# -----------
# Parameters:
# -----------

# these args are just used to put re-runs of the analysis in different output directories 
ANALYSIS_NAME="no_expr_pc_covars"
RESOLUTION="major_cell_types"

i=${SGE_TASK_ID};
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)

# --------
# Get significant SNPs from GeNA output
# --------
# Reference: https://github.com/immunogenomics/GeNA-applied/blob/main/run_gwas/02.munge_geno_sigsnps.sh 

# define inputs
# GWAS=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}/${CELLTYPE}/${ANALYSIS_NAME}/GeNA_sumstats.txt
SNP_LIST="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}/${CELLTYPE}/${ANALYSIS_NAME}/GeNA_snp_list_sig_MAF_0.05.tsv"
GENOTYPES=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15
ALT_ALLELES="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}/${CELLTYPE}/${ANALYSIS_NAME}/GeNA_snp_list_sig_MAF_0.05_alt_alleles.tsv"

# define outputs 
OUT_PLINK=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}/${CELLTYPE}/${ANALYSIS_NAME}/GeNA_sig_snps

# 1. export plinks to plink.raw tsv 

plink2 --pfile ${GENOTYPES} \
    --extract ${SNP_LIST} \
    --threads 10 \
    --export A \
    --export-allele ${ALT_ALLELES} \
    --out ${OUT_PLINK}

# 2. Run the python pruning script 

python /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA/08_prune_sig_snps.py ${CELLTYPE} ${RESOLUTION} ${ANALYSIS_NAME}


############ Graveyard 🪦 ##############

# permformed by 06_GeNA_plotting.R
# p_thresh=5e-08

# # 2. Generate a plink file containing only the significant SNPs

# plink2 --pfile ${GENOTYPES} \
#     --extract ${SNP_LIST} \
#     --threads 10 \
#     --make-pgen \
#     --out ${OUT_PLINK}

# # 3. generate a VCF file containing only the significant SNPs 

# plink2 --pfile ${GENOTYPES} \
#     --extract ${SNP_LIST} \
#     --threads 10 \
#     --export vcf bgz \
#     --out ${OUT_PLINK}


# AT the moment, just using this one for downstream analysis. 
# 4. Genotype file for all significant SNPS 
