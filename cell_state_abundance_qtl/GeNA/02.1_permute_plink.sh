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
