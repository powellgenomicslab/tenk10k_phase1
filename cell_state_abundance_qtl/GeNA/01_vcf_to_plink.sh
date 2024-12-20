#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=10G
#$ -pe smp 10
#$ -N vcf_to_plink
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/vcf_to_plink.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/vcf_to_plink.stdout
#$ -m ae
#$ -t 1-22
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate mastectomy-env

cd ${TMPDIR}

i=${SGE_TASK_ID};

IN_VCF=$(ls /directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/genotypes/december2024_freeze/chr*_common_variants.vcf.bgz | awk -v task_id=${SGE_TASK_ID} 'NR==task_id')

FILE_BASENAME=$(basename ${IN_VCF})
PFILE="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/${FILE_BASENAME%.vcf.bgz}_standard_chr"
TMP_VCF=${FILE_BASENAME}_standard_chr.vcf.gz

# keep only standard chromosomes
echo "bcftools view --regions chr{1..22},chrX,chrY,chrM ${IN_VCF} -Oz -o ${TMP_VCF}"
bcftools view --regions $(printf "chr%s," {1..22})chrX,chrY,chrM ${IN_VCF} -Oz -o ${TMP_VCF}

# keep only SNPs 
echo "plink2 --threads 10 --vcf ${TMP_VCF} --make-pgen --allow-extra-chr --snps-only --out ${PFILE}_snps"
plink2 --threads 10 --vcf ${TMP_VCF} --make-pgen --allow-extra-chr --snps-only --out ${PFILE}_snps
