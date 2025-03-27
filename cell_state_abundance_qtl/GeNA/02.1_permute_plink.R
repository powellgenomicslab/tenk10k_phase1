# mm activate mastectomy-env


# Generate permuted / shuffled cpg ids
# making sure the permuted sample ids are 
library(tidyverse)

# sample_meta <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/saige-qtl_tenk10k-genome-2-3-eur_input_files_241210_covariates_sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv")

# get sample id's present in the vcf file
system("bcftools query -l /directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/genotypes/december2024_freeze/chr1_common_variants.vcf.bgz > /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/sample_ids_vcf.txt")
sample_ids_vcf <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/sample_ids_vcf.txt", col_names = "sample_id")

# metadata <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/298_libraries_cell_metadata_filtered_min1000genes.csv")

# sample_ids_single_cell <- metadata %>%
#     select(cpg_id) %>%
#     distinct()

set.seed(123)

permutation_mapping <- sample_ids_vcf %>%
    select(sample_id) %>%
    mutate(
        sample_id_perm1 = sample(sample_id),
        sample_id_perm2 = sample(sample_id),
        sample_id_perm3 = sample(sample_id),
    )

write_csv(permutation_mapping, "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/sample_ids_vcf_with_perm_ids.txt")



# table(duplicated(permutation_mapping$sample_id_perm2))





# Analysis to check how the permuted ID's line up with ids in the single cell object and in the genotype VCF


# # 1 ID in VCF but not in sample_meta
# setdiff(sample_ids_vcf, sample_meta$sample_id)

# # all single cell samples in sample_ids_vcf 
# setdiff(sample_ids_single_cell$cpg_id, sample_ids_vcf)

# # get the permutation mappings that correspond to the single cell real ids 
# sample_meta_in_cs <- sample_meta %>%
#     filter(sample_id %in% sample_ids_single_cell$cpg_id)

# # 309 ID's in perm0 are not included in the genotype VCF file
# setdiff(sample_meta_in_cs$sample_perm0, sample_ids_vcf) %>%
#     length()

# NOTE:
# need to redo the permutation so that *all* permuted ID's are actually in the genotype VCF 
# currently these might not be getting used in the permuted analysis
# possible that issues with permutation are due to ~300ish samples with no genotype info 
