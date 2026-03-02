# mm activate mastectomy-env

# Generate permuted / shuffled cpg ids
# making sure the permuted sample ids are
library(tidyverse)
library(glue)
# sample_meta <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/saige-qtl_tenk10k-genome-2-3-eur_input_files_241210_covariates_sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv")

# get sample id's present in the vcf file
system("bcftools query -l /directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/genotypes/december2024_freeze/chr1_common_variants.vcf.bgz > /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/sample_ids_vcf.txt")
sample_ids_vcf <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/sample_ids_vcf.txt", col_names = "sample_id")

set.seed(123)
# used perm1 for the actual analysis // calibration
permutation_mapping <- sample_ids_vcf %>%
    select(sample_id) %>%
    mutate(
        sample_id_perm1 = sample(sample_id),
        sample_id_perm2 = sample(sample_id),
        sample_id_perm3 = sample(sample_id),
    )

write_csv(permutation_mapping, "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/sample_ids_vcf_with_perm_ids.txt")

