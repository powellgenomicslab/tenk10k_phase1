# PURPOSE: Create a master metadata file to keep track of sample naming, pooling and other important stuff
library(tidyverse)
library(janitor)
library(RSQLite)

# read cpg mappings and donor pools info
cpg_map_file <- "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/str_sample-sex-mapping_sample_karyotype_sex_mapping.csv"
donor_pools_bioheart_file <- "/share/ScratchGeneral/anncuo/tenk10k/donor_info/TenK10K_BioHeart_pool_info_v3.csv"
donor_pools_tob_file <- "/share/ScratchGeneral/anncuo/tenk10k/donor_info/TOB_TenK10K_pool_info_v3.csv"
onek1k_file_name <- "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/OneK1K/published_data/Database/OneK1K_MetadataDB.db"

# CPG ID's (from WGS)
cpg_mapping <- read_csv(cpg_map_file) %>%
    rename(cpg_id = s)

# Bioheart
donor_pools_bioheart <- read_csv(donor_pools_bioheart_file) %>%
    janitor::clean_names() %>%
    mutate(bioheart_cohort = TRUE) %>%
    rename(
        bioheart_id = sample_id,
        bioheart_pool = pool
    )

# TOB (from onek1k)
donor_pools_tob <- read_csv(donor_pools_tob_file) %>%
    janitor::clean_names() %>%
    mutate(external_id_modified = tob_id %>% str_remove("-0")) %>% # clean tob mapping id so the data can be joined
    mutate(tob_cohort = TRUE) %>%
    mutate(onek1k_cohort = TRUE)

# Onek1k metadata
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver, dbname = onek1k_file_name)
onek1k_metadata <- dbReadTable(db, "RECRUITED_DONORS") %>%
    tibble() %>%
    janitor::clean_names()

nrow(onek1k_metadata)
nrow(donor_pools_tob)
colnames(donor_pools_tob)

# join the TOB data
all_tob <- onek1k_metadata %>%
    full_join(donor_pools_tob, by = "tob_id")

table(duplicated(donor_pools_tob$tob_id))

master_table <- cpg_mapping %>%
    mutate(external_id_modified = external_id %>% str_remove("-PBMC")) %>%
    full_join(donor_pools_tob, by = "external_id_modified") %>%
    full_join(donor_pools_bioheart, by = c("external_id_modified" = "bioheart_id")) %>%
    mutate(tenk10k_pool = coalesce(tenk10k_pool, bioheart_pool))

nrow(donor_pools_tob)
sum(!is.na(master_table$tob_id))

master_table %>%
    write_csv("data_processing/tenk10k_donor_pool_info.csv")

# TODO: add which cohort the donor is from (TOB, bioheart, onek1k)
# TODO: add which pool each cohort is from
# TODO: add sequencing batches - each row will correspond to the sequencing batch i.e. will include the repeats
cpg_mapping %>% write_csv("data_processing/str_sample-sex-mapping_sample_karyotype_sex_mapping.csv")

onek1k_metadata %>% write_csv("data_processing/OneK1K_Metadata_donors.csv")
