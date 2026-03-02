#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=50G
#$ -l tmp_requested=40G
#$ -N coloc_prep_geno
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc_prep_geno.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/coloc_prep_geno.stdout
#$ -t 1-8
#$ -m ae
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate mastectomy-env

ANALYSIS_NAME="no_expr_pc_covars"
RESOLUTION="major_cell_types"

i=${SGE_TASK_ID};
CELLTYPE=$(sed "${i}q;d" /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/major_cell_types.txt)

GENOTYPES_ALL=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15
SNP_DIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}/${CELLTYPE}/${ANALYSIS_NAME}/lead_snps_ciswindow/

# 1. for each lead csaQTL snp, create a genotype plink file containing the snps within its cis-window 

LEAD_SNPS_CHR_POS="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}/${CELLTYPE}/${ANALYSIS_NAME}/GeNA_sumstats_lead_snps_MAF_0.05.snps"

# make an array containing all of the lead SNPs
lead_snps=()
while IFS= read -r line; do
    lead_snps+=("$line")
done < $LEAD_SNPS_CHR_POS

# iterate through the list of lead SNPs for the cell type,

for lead_snp in "${lead_snps[@]}"; do
    echo "Processing $lead_snp in $CELLTYPE"
    CIS_SNPS_LIST="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}/${CELLTYPE}/${ANALYSIS_NAME}/lead_snps_ciswindow/GeNA_sumstats_${lead_snp}_cis_snps.snps"
    GENOTYPES_CIS_SNPS=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/${RESOLUTION}/${CELLTYPE}/${ANALYSIS_NAME}/lead_snps_ciswindow/${lead_snp}_cis_snps
    # create plink file containing genotypes for the SNPs in cis window around the lead snp 
    plink2 --pfile ${GENOTYPES_ALL} \
        --extract ${CIS_SNPS_LIST} \
        --make-pgen \
        --out ${GENOTYPES_CIS_SNPS}

    # for each lead SNP, run PLINK to generate fixed-phenotype summary statistics which are compatible with coloc analysis 
    # Refer to https://github.com/immunogenomics/GeNA-applied/blob/main/run_gwas/06.local_test_fixed_pheno.sh
    TRAIT_FILE="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/${RESOLUTION}/${ANALYSIS_NAME}/${CELLTYPE}_sample_pheno_for_plink.tsv"
    mkdir -p /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/plink_fixed_pheno/${CELLTYPE}/
    OUTPUT_FIXED_PHENOTPE=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/coloc/plink_fixed_pheno/${CELLTYPE}/gwas_local_${lead_snp}_cis
    plink2 --pfile ${GENOTYPES_CIS_SNPS} \
        --pheno ${TRAIT_FILE} \
        --glm "allow-no-covars" \
        --prune \
        --pheno-name spheno_${lead_snp} \
        --out ${OUTPUT_FIXED_PHENOTPE}

done

