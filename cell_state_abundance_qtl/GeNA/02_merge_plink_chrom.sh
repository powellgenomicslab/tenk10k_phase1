#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -r yes
#$ -l mem_requested=10G
#$ -pe smp 10
#$ -N merge_plink
#$ -q short.q
#$ -e /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/merge_plink.stderr
#$ -o /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/logs/merge_plink.stdout
#$ -m ae
#$ -M b.bowen@garvan.org.au

. /home/${USER}/micromamba/etc/profile.d/micromamba.sh
micromamba activate mastectomy-env

CHR1_BFILE=/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/chr1_common_variants_standard_chr_snps

# merge indiv chromosome plink files into one big file 
plink2 --pfile /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/chr1_common_variants_standard_chr_snps \
    --pmerge-list /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/plink_merge_chr.txt \
    --threads 10 \
    --out /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_snps\

