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
PERM_MAPPING=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/sample_ids_vcf_with_perm_ids.txt
PFILE_MERGED="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15"

# 1. create mapping between real and permuted ID's 
awk -F, 'NR > 1 {OFS="\t"; print $1, $2}' ${PERM_MAPPING} > ${OUTDIR}/perm1_map.txt

# 2. Reheader the plink file using the permuted ID's 
plink2 --pfile ${PFILE_MERGED} \
    --update-ids ${OUTDIR}/perm1_map.txt \
    --make-pgen \
    --out ${OUTDIR}/merged_common_variants_standard_chr_geno_0.15_singlecell_cohort_perm1 \
    --threads 4


######################################################
# Permutation testing for single variant
######################################################

# # Generate a bunch of perumutations for just 1 variant of interest 

# VARIANT="15:39687137:C:T"
# OUTDIR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted/variant_permutations
# PFILE_MERGED="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15"
# GENO_VAR=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_single_variant/${VARIANT}_standard_chr_geno_0.15
# PERM_MAPPING=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/sample_ids_vcf_with_perm_ids_1000.txt
# VAR_FILE=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted/variant_permutations/variant_list.txt

# #  create txt file containing the variant name 
# echo ${VARIANT} > ${VAR_FILE}
# # create plink file containing just the variant of interest 
# plink2 --pfile ${PFILE_MERGED} \
#     --extract ${VAR_FILE} \
#     --make-pgen \
#     --out ${GENO_VAR} \
#     --threads 4

# for i in {1..1000}; do 
#     echo " Creating permuted data for permutation # ${i}"
#     # create a txt file with the mapping for each permutation 
#     awk -F,  -v col="$i" 'NR > 1 {OFS="\t"; print $1, $(col + 1)}' ${PERM_MAPPING} > ${OUTDIR}/perm${i}_map.txt
#     # shuffle the sample ids in the plink file according to that mapping 
    
#     # 2. Reheader the plink file using the permuted ID's 
#     plink2 --pfile ${GENO_VAR} \
#         --update-ids ${OUTDIR}/perm${i}_map.txt \
#         --make-pgen \
#         --out ${OUTDIR}/${VARIANT}_perm${i} \
#         --threads 4
# done