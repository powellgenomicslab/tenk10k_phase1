#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=10G
#$ -l tmp_requested=200G
#$ -pe smp 4
#$ -N permute_plink
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/permute_plink.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/permute_plink.stdout
#$ -m ae
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate mastectomy-env

OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted

SAMPLE_META="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/saige-qtl_tenk10k-genome-2-3-eur_input_files_241210_covariates_sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv"

PFILE_MERGED="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15"

# 1. create list of samples to keep, these are the samples in both the single cell and wgs cohort
# awk -F, 'NR > 1 {print $1}' ${SAMPLE_META} > ${OUTDIR}/samples_keep.txt

# 2. create mapping between real and permuted ID's 
awk -F, 'NR > 1 {OFS="\t"; print $1, $20}' ${SAMPLE_META} > ${OUTDIR}/perm0_map.txt

# 3. Filter plink to only include the sc cohort samples 
# plink2 --pfile ${PFILE_MERGED}\
#     --keep ${OUTDIR}/samples_keep.txt \
#     --make-pgen \
#     --out ${TMPDIR}/merged_common_variants_standard_chr_geno_0.15_subset \
#     --threads 4

# PFILE_SUBSET=${TMPDIR}/merged_common_variants_standard_chr_geno_0.15_subset


# 4. Reheader the plink file using the permuted ID's 
plink2 --pfile ${PFILE_MERGED} \
    --update-ids ${OUTDIR}/perm0_map.txt \
    --make-pgen \
    --out ${OUTDIR}/merged_common_variants_standard_chr_geno_0.15_singlecell_cohort_perm0 \
    --threads 4
