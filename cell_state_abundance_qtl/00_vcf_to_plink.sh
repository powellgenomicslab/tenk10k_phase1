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
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate gena-env

cd ${TMPDIR}

IN_VCF=/share/ScratchGeneral/anncuo/tenk10k/data_processing/genotypes/Merged_MAF0.05_hg38_chr.vcf.gz
FILE_BASENAME=$(basename ${IN_VCF})
PFILE="/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/${FILE_BASENAME%.vcf.gz}"
TMP_VCF=Merged_MAF0.05_R2_0.3_hg38_standard_chr.vcf.gz

# keep only standard chromosomes
echo "bcftools view --regions chr{1..22},chrX,chrY,chrM ${IN_VCF} -Oz -o ${TMP_VCF}"
bcftools view --regions $(printf "chr%s," {1..22})chrX,chrY,chrM ${IN_VCF} -Oz -o ${TMP_VCF}

# keep only SNPs 
echo "plink2 --threads 10 --vcf ${TMP_VCF} --make-pgen --allow-extra-chr --snps-only --out ${PFILE}_snps"
plink2 --threads 10 --vcf ${TMP_VCF} --make-pgen --allow-extra-chr --snps-only --out ${PFILE}_snps
