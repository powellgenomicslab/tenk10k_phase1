# PURPOSE: Create a master metadata file to keep track of sample naming, pooling and other important stuff
library(tidyverse)
library(janitor)
library(RSQLite)
library(glue)

# TODO: figure out why there are only 211 pools; should be 300
# I think this is because the metadata has the pooling info but many pools were repeated / split across flowcells etc and the 211 pools doesn't take that into account whereas th 300 does

# read cpg mappings and donor pools info
cpg_map_file <- "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/str_sample-sex-mapping_sample_karyotype_sex_mapping.csv"
cpg_map_file_updated_bh <- "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/metadata/bioheart_ct_ids_to_cpg_ids_20241113_n1419.csv"
donor_pools_bioheart_file <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/TenK10K_BioHeart_pool_info_v5.csv"
donor_pools_tob_file <- "/share/ScratchGeneral/anncuo/tenk10k/donor_info/TOB_TenK10K_pool_info_v3.csv"
onek1k_file_name <- "/directflow/SCCGGroupShare/projects/data/experimental_data/projects/OneK1K/published_data/Database/OneK1K_MetadataDB.db"

# CPG ID's (from WGS)
cpg_mapping <- read_csv(cpg_map_file) %>%
    rename(cpg_id = s) %>%
    select(-sex_karyotype) %>%
    mutate(notes = if_else(str_detect(external_id, "-PBMC"), glue("Alternate external ID: {external_id}. -PBMC suffix indicates that these samples had WGS re-run due to low quality."), NA)) %>%
    mutate(external_id = str_remove(external_id, "-PBMC"))

cpg_mapping_bioheart <- read_csv(cpg_map_file_updated_bh) %>%
    select(-1) %>%
    rename("external_id" = CT_ID, "cpg_id" = CPG_ID) %>%
    mutate(notes = if_else(str_detect(external_id, "-PBMC"), glue("Alternate external ID: {external_id}. -PBMC suffix indicates that these samples had WGS re-run due to low quality."), NA)) %>%
    mutate(external_id = str_remove(external_id, "-PBMC"))

cpg_mapping_complete <- bind_rows(cpg_mapping, cpg_mapping_bioheart) %>%
    distinct()

# NOTE: confirmed no duplicates
# cpg_mapping_complete %>%
#     pull(cpg_id) %>%
#     duplicated() %>%
#     table()

# cpg_mapping_bioheart %>% filter(external_id == "CT_2924")

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
# NOTE: only 38 pools - is this correct??
donor_pools_tob <- read_csv(donor_pools_tob_file) %>%
    janitor::clean_names() %>%
    # mutate(external_id_modified = tob_id %>% str_remove("-0")) %>% # clean tob mapping id so the data can be joined
    mutate(tob_cohort = TRUE) %>%
    mutate(onek1k_cohort = TRUE)

# donor_pools_tob %>% filter(tob_id %in% c("TOB-01003", "TOB-01002"))

# Onek1k metadata
sqlite.driver <- dbDriver("SQLite")
db <- dbConnect(sqlite.driver, dbname = onek1k_file_name)
tob_id_mapping <- dbReadTable(db, "RECRUITED_DONORS") %>%
    tibble() %>%
    janitor::clean_names() %>%
    select(patient_id, tob_id, notes, iid, study_id) %>%
    mutate(
        notes = if_else(notes == "NA", NA, notes),
        tob_id = if_else(tob_id == "NA", NA, tob_id)
    ) %>%
    distinct() %>%
    drop_na(tob_id)

# tob_id_mapping %>% filter(tob_id %in% c("TOB-01003", "TOB-01002"))

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

# donor_pools_tob %>% filter(tob_id %in% c("TOB-01003", "TOB-01002"))

# donor_pools_tob %>%
#     filter(is.na(patient_id)) %>%
#     print(n = 100)

# donor_pools_tob[donor_pools_tob$tob_id %in% donor_pools_tob$tob_id[duplicated(donor_pools_tob$tob_id)], ] %>% arrange(tob_id, tenk10k_pool)

# NOTES:
# 1. Full join BH pool info with TOB/CPG mappings: includes BH samples even if they don't currently have CPG-ids
# 2: there is one duplicated TOB sample that is not marked as a repeat, although it maps to 2 different pools so probably is just an un-marked repeat sample
# 3. keeps the notes in there as they are useful for keeping track of the repeated samples

master_table <- cpg_mapping_complete %>%
    full_join(donor_pools_tob, by = c("external_id" = "patient_id")) %>% # full join so that we retain TOB libraries without CPG id's
    mutate(notes = coalesce(notes.x, notes.y)) %>%
    select(-notes.x, -notes.y) %>%
    # # full join so that samples without CPG_ids are included (while waiting for WGS)
    full_join(donor_pools_bioheart, by = c("external_id" = "bioheart_id")) %>%
    mutate(notes = coalesce(notes.x, notes.y)) %>%
    select(-notes.x, -notes.y) %>%
    mutate(tenk10k_pool = coalesce(tenk10k_pool, bioheart_pool)) %>%
    select(cpg_id, external_id, tob_id, study_id, tenk10k_pool, tob_cohort, bioheart_cohort, notes) %>%
    group_by(cpg_id, external_id, tob_id, study_id, tenk10k_pool, tob_cohort, bioheart_cohort) %>%
    filter(!(n() > 1 & is.na(notes))) %>% # remove the duplicate ID's where NOTES field is missing
    ungroup()

master_table %>%
    write_csv("data_processing/tenk10k_donor_pool_info.csv")


# check that new sample CPG-ID mappings are available

new_cpg_ids <- c(
    # 240223
    "S0147",
    "S0154",
    "S0155",
    "S0156",
    "S0157",
    "S0158",
    "S0159",
    "S0160",
    "S0161",
    "S0162",
    "S0163",
    "S0164",
    "S0165",
    "S0167",
    "S0168",
    "S0169",
    "S0170",
    "S0017",
    "S0019",
    "S0166",
    "S0171",
    "S0172",
    "S0173",
    "S0174",
    "S0175",
    "S0176",
    "S0178",
    "S0179",
    "S0180",
    "S0181",
    "S0182",
    "S0183",
    "S0184",
    # 240501
    "S0017a",
    "S0019a",
    "S0166",
    "S0171",
    "S0172",
    "S0173",
    "S0174",
    "S0175",
    "S0176",
    "S0178",
    "S0179",
    "S0180",
    "S0181",
    "S0182",
    "S0183",
    "S0184",
    # 240524
    "S0177",
    "S0185",
    "S0186",
    "S0187",
    "S0188_re1",
    "S0189",
    "S0190",
    "S0191",
    "S0192",
    "S0193",
    "S0194",
    "S0195",
    "S0196",
    "S0197",
    "S0198",
    "S0199",
    "S0200",
    "S0201",
    "S0202",
    "S0203",
    "S0204",
    "S0205",
    "S0206",
    "S0207",
    "S0208",
    "S0209",
    "S0210"
)

# Missing NEW CPG id's confirmed to be QC fails etc.
#  likely to be topped up after phase 1?
new_cpg_ids[!new_cpg_ids %in% master_table$tenk10k_pool]
missing_cpg <- master_table %>%
    filter(tenk10k_pool %in% new_cpg_ids) %>%
    filter(is.na(cpg_id)) %>%
    pull(external_id)
# arrange(tenk10k_pool) %>%
# print(n = 100)

missing_cpg

# check ID duplication ----
duplicated_ids <- master_table %>%
    filter(duplicated(external_id)) %>%
    select(external_id) %>%
    distinct() %>%
    pull(external_id)
master_table %>%
    filter(external_id %in% duplicated_ids) %>%
    print(n = 100)

# # check ID duplication:
# duplicated_ids_tob_id <- master_table %>%
#     filter(duplicated(tob_id)) %>%
#     select(tob_id) %>%
#     drop_na() %>%
#     distinct() %>%
#     pull(tob_id)
# master_table %>%
#     filter(tob_id %in% duplicated_ids_tob_id) %>%
#     print(n = 100)

# master_table %>% filter(external_id %in% c("TOB1003", "TOB1002"))


cpg_map_file <- read_csv("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/metadata/bioheart_ct_ids_to_cpg_ids_20241113_n1419.csv")

vcf_individuals_subset <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cpg_ids_bioheart_367new_oct2024_demux_vcf_common_variants.txt")
vcf_individuals <- read_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cpg_ids_tob_bioheart_combined_nov2024_demux_vcf_common_variants.txt")

missing_ids <- cpg_map_file$CPG_ID[!cpg_map_file$CPG_ID %in% vcf_individuals]
# missing_ids2 <- cpg_map_file$CPG_ID[!cpg_map_file$CPG_ID %in% vcf_individuals_subset]

cpg_map_file$CPG_ID[!cpg_map_file$CPG_ID %in% vcf_individuals] %>% write_lines("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/cpg_ids_missing_from_vcf.txt")

table(missing_ids %in% vcf_individuals_subset)
table(missing_ids %in% vcf_individuals)
