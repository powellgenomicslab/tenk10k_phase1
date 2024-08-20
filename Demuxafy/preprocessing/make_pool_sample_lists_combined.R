library(tidyverse)
library(data.table)
library(glue)

# USAGE: creates a master metadata sheet for both * BioHEART and TOB * samples.
# Includes donor information (which BioHEART IDs are in which pool)

# NOTE:
# The input file "tenk10k_donor_pool_info.csv" that this script relies on can be generated with create_master_metadata.R
# May need to update the input files to get the latest updates in the google sheets
# /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/Demuxafy/preprocessing/create_master_metadata.R

# files in this directory have the list of samples within each pool in tenk10k
pools_dir <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/samples_in_pools_combined/"

## get donors in pools for both bioheart and TOB
# if this shows 0 donors in a pool check that the TenK10K_BioHeart_pool_info table has been updated, and check that the column names have not been changed
get_n_donors_all_pools <- function(pool_list_suffix, pool_list_path, save = FALSE, pools_dir) {
    # get list of all samples in the sequencing library
    all_samples_in_pool <- fread(glue("{pool_list_path}/cellranger_outs_{pool_list_suffix}.txt"), header = FALSE)

    # creates a file for each sequencing library showing the donors in that sample pool
    get_donors_in_libraries <- function(pool, pools_dir, save) {
        # strip text from pool name indicated "repeat" status
        tenk10k_pool <- str_remove(pool, pattern = "_re$|_repeat$|a$|b$|c$|d$|_re1$")
        print(paste0("Input pool ID: ", pool, "  ", "tenk pool ID: ", tenk10k_pool))

        # read in the master pool info sheet
        tenk10k_all_pool_info <- fread("data_processing/tenk10k_donor_pool_info.csv") %>%
            select(external_id, cpg_id, tenk10k_pool)

        donors_in_pool <- tenk10k_all_pool_info %>%
            filter(tenk10k_pool == {{ tenk10k_pool }}) %>%
            select(external_id)

        colnames(donors_in_pool) <- c()
        print(donors_in_pool)

        new_file <- paste0(pools_dir, pool, ".tsv") # saves with the input pool name - NOT the tenk10k pool name
        if (save) {
            fwrite(donors_in_pool, new_file, sep = "\t") # create a file with the correct sequencing library name
        } else {
            return(donors_in_pool)
        }
    }

    # this function gets the number of donors for a given pool
    get_n_donors_pool <- function(pool) {
        donors_in_pool_file <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/samples_in_pools_combined/{pool}.tsv") # file generated in the above function

        # get number of rows (donors) that are in each sequencing library
        donors_in_pool <- fread(donors_in_pool_file, header = FALSE) %>%
            nrow()
        return(data.frame(pool = pool, donors_in_pool = donors_in_pool))
    }

    # for each sample/sequencing library/pool:
    output <- all_samples_in_pool %>%
        pull(V1) %>% # 1. creates file containing the donors in each sequencig library
        purrr::walk(\(x) get_donors_in_libraries(
            pool = x,
            pools_dir = pools_dir,
            save = save
        )) %>% # 2. counts how many donors there are and saves a table with total n donors for each of the sequencing libraries
        purrr::map(
            \(x) get_n_donors_pool(pool = x)
        ) %>%
        list_rbind() # make a df with donor numbers for all the samples contained in the cellranger outs

    if (save) {
        # print("saving...")
        output_path <- glue("{pool_list_path}/libraries_nsamples_{pool_list_suffix}.txt")
        write_tsv(output, output_path, col_names = FALSE)
    }
    return(output)
}


# PREVIOUS RUNS:
# get_n_donors_all_pools(
#     pool_list_suffix = "240501",
#     pool_list_path = "data_processing",
#     save = TRUE,
#     pools_dir = pools_dir
# )
# get_n_donors_all_pools(
#     pool_list_suffix = "240524",
#     pool_list_path = "data_processing",
#     save = TRUE,
#     pools_dir = pools_dir
# )

get_n_donors_all_pools(
    pool_list_suffix = "231213",
    pool_list_path = "data_processing",
    save = TRUE,
    pools_dir = pools_dir
)
get_n_donors_all_pools(
    pool_list_suffix = "231214",
    pool_list_path = "data_processing",
    save = TRUE,
    pools_dir = pools_dir
)
# NOTE: ^ script does not account for the accidental maxi-pools S0034-37a and S0034-37b, need to manually create these?






#### Run this function to get all the samples in pools for the older batches ----
tenk10k_all_pool_info <- fread("data_processing/tenk10k_donor_pool_info.csv") %>%
    select(external_id, cpg_id, tenk10k_pool)

# creates a file for each sequencing library showing the donors in that sample pool
get_donors_in_libraries <- function(pool, pools_dir, save) {
    # strip text from pool name indicated "repeat" status
    tenk10k_pool <- str_remove(pool, pattern = "_re$|_repeat$|a$|b$|c$|d$|_re1$")
    print(paste0("Input pool ID: ", pool, "  ", "tenk pool ID: ", tenk10k_pool))

    # read in the master pool info sheet
    tenk10k_all_pool_info <- fread("data_processing/tenk10k_donor_pool_info.csv") %>%
        select(external_id, cpg_id, tenk10k_pool)

    donors_in_pool <- tenk10k_all_pool_info %>%
        filter(tenk10k_pool == {{ tenk10k_pool }}) %>%
        select(external_id)

    colnames(donors_in_pool) <- c()
    print(donors_in_pool)

    new_file <- paste0(pools_dir, pool, ".tsv") # saves with the input pool name - NOT the tenk10k pool name
    if (save) {
        fwrite(donors_in_pool, new_file, sep = "\t") # create a file with the correct sequencing library name
    } else {
        return(donors_in_pool)
    }
}

all_pools <- unique(tenk10k_all_pool_info$tenk10k_pool)

all_pools %>% map(\(pool) get_donors_in_libraries(pool = pool, pools_dir = pools_dir, save = TRUE))

##### -----

# make a list of all the CPG ID's for vireo all samples run

all_cpg <- tenk10k_all_pool_info %>%
    select(cpg_id) %>%
    filter(!is.na(cpg_id))
fwrite(all_cpg, "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/all_cpg_ids.txt", col.names = FALSE)

tob_all_pool_info <- fread("data_processing/tenk10k_donor_pool_info.csv") %>% filter(tob_cohort == TRUE)

all_cpg_tob <- tob_all_pool_info %>%
    select(cpg_id) %>%
    filter(!is.na(cpg_id))
fwrite(all_cpg_tob, "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/all_cpg_ids_tob.txt", col.names = FALSE)
