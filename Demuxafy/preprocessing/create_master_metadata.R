# PURPOSE: Create a master metadata file to keep track of sample naming, pooling and other important stuff
library(tidyverse)
library(janitor)
library(RSQLite)
library(glue)

# read cpg mappings and donor pools info
cpg_map_file <- "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/str_sample-sex-mapping_sample_karyotype_sex_mapping.csv"
donor_pools_bioheart_file <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/TenK10K_BioHeart_pool_info_v5.csv"
donor_pools_tob_file <- "/share/ScratchGeneral/anncuo/tenk10k/donor_info/TOB_TenK10K_pool_info_v3.csv"
onek1k_file_name <- "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/OneK1K/published_data/Database/OneK1K_MetadataDB.db"

# CPG ID's (from WGS)
cpg_mapping <- read_csv(cpg_map_file) %>%
    rename(cpg_id = s) %>%
    select(-sex_karyotype) %>%
    mutate(notes = if_else(str_detect(external_id, "-PBMC"), glue("Alternate external ID: {external_id}. -PBMC suffix indicates that these samples had WGS re-run due to low quality."), NA)) %>%
    mutate(external_id = str_remove(external_id, "-PBMC"))

# Bioheart
donor_pools_bioheart <- read_csv(donor_pools_bioheart_file) %>%
    janitor::clean_names() %>%
    mutate(bioheart_cohort = TRUE) %>%
    rename(
        bioheart_id = bio_heart_id,
        bioheart_pool = pool
    ) %>%
    select(-cpg_id)

# TOB (from onek1k)
donor_pools_tob <- read_csv(donor_pools_tob_file) %>%
    janitor::clean_names() %>%
    # mutate(external_id_modified = tob_id %>% str_remove("-0")) %>% # clean tob mapping id so the data can be joined
    mutate(tob_cohort = TRUE) %>%
    mutate(onek1k_cohort = TRUE)

donor_pools_tob %>% filter(tob_id %in% c("TOB-01003", "TOB-01002"))

# Onek1k metadata
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver, dbname = onek1k_file_name)
tob_id_mapping <- dbReadTable(db, "RECRUITED_DONORS") %>%
    tibble() %>%
    janitor::clean_names() %>%
    select(patient_id, tob_id, notes, iid) %>%
    mutate(
        notes = if_else(notes == "NA", NA, notes),
        tob_id = if_else(tob_id == "NA", NA, tob_id)
    ) %>%
    distinct() %>%
    drop_na(tob_id)

tob_id_mapping %>% filter(tob_id %in% c("TOB-01003", "TOB-01002"))

# tob_id_mapping[tob_id_mapping$tob_id %in% tob_id_mapping$tob_id[duplicated(tob_id_mapping$tob_id)], ] %>% arrange(patient_id)

# NOTE: FID column - there are a lot of samples that have 2 FIDs ~ 16 donors

# tob_id_mapping_old <- dbReadTable(db, "RECRUITED_DONORS") %>%
#     tibble() %>%
#     janitor::clean_names() %>%
#     select(patient_id, tob_id) %>%
#     distinct()

# join the TOB/onek1k data
donor_pools_tob <- donor_pools_tob %>%
    left_join(tob_id_mapping, by = c("tob_id")) %>%
    mutate(notes = coalesce(notes.x, notes.y)) %>%
    select(-notes.x, -notes.y) %>%
    group_by(tob_id, tenk10k_pool) %>%
    filter(!(n() > 1 & is.na(notes))) %>% # remove the duplicate ID's where NOTES field is missing
    ungroup()

donor_pools_tob %>% filter(tob_id %in% c("TOB-01003", "TOB-01002"))

# donor_pools_tob %>%
#     filter(is.na(patient_id)) %>%
#     print(n = 100)

# donor_pools_tob[donor_pools_tob$tob_id %in% donor_pools_tob$tob_id[duplicated(donor_pools_tob$tob_id)], ] %>% arrange(tob_id, tenk10k_pool)

# NOTES:
# 1. Full join BH pool info with TOB/CPG mappings: includes BH samples even if they don't currently have CPG-ids
# 2: there is one duplicated TOB sample that is not marked as a repeat, although it maps to 2 different pools so probably is just an un-marked repeat sample
# 3. keeps the notes in there as they are useful for keeping track of the repeated samples

master_table <- cpg_mapping %>%
    left_join(donor_pools_tob, by = c("external_id" = "patient_id")) %>%
    mutate(notes = coalesce(notes.x, notes.y)) %>%
    select(-notes.x, -notes.y) %>%
    # # full join so that samples without CPG_ids are included (while waiting for WGS)
    full_join(donor_pools_bioheart, by = c("external_id" = "bioheart_id")) %>%
    mutate(notes = coalesce(notes.x, notes.y)) %>%
    select(-notes.x, -notes.y) %>%
    mutate(tenk10k_pool = coalesce(tenk10k_pool, bioheart_pool)) %>%
    select(cpg_id, external_id, tob_id, tenk10k_pool, tob_cohort, bioheart_cohort, notes) %>%
    group_by(cpg_id, external_id, tob_id, tenk10k_pool, tob_cohort, bioheart_cohort) %>%
    filter(!(n() > 1 & is.na(notes))) %>% # remove the duplicate ID's where NOTES field is missing 
    ungroup()

master_table %>%
    write_csv("data_processing/tenk10k_donor_pool_info.csv")

# check ID duplication:
duplicated_ids <- master_table %>%
    filter(duplicated(external_id)) %>%
    select(external_id) %>%
    distinct() %>%
    pull(external_id)
master_table %>%
    filter(external_id %in% duplicated_ids) %>%
    print(n = 100)

# check ID duplication:
duplicated_ids_tob_id <- master_table %>%
    filter(duplicated(tob_id)) %>%
    select(tob_id) %>%
    drop_na() %>% 
    distinct() %>%
    pull(tob_id)
master_table %>%
    filter(tob_id %in% duplicated_ids_tob_id) %>%
    print(n = 100)

master_table %>% filter(external_id %in% c("TOB1003", "TOB1002"))
