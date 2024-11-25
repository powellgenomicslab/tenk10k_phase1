library(tidyverse)
library(data.table)
library(glue)

# USAGE: creates a master metadata sheet for both * BioHEART and TOB * samples.
# Includes donor information (which BioHEART IDs are in which pool)

# NOTE:
# The input file "tenk10k_donor_pool_info.csv" that this script relies on can be generated with create_master_metadata.R
# May need to update the input files to get the latest updates in the google sheets
# /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/Demuxafy/preprocessing/create_master_metadata.R

################################################################
################### Most BIOheart + TOB ########################
################################################################

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
        
        # for debugging 
        tenk10k_all_pool_info %>%
            filter(tenk10k_pool == {{ tenk10k_pool }}) %>%
            print()
        
        # get all cpg ids in the pool (may not be cpg ids for all scRNA-seq donors)
        donors_in_pool <- tenk10k_all_pool_info %>%
            filter(tenk10k_pool == {{ tenk10k_pool }}) %>%
            select(cpg_id) %>%
            filter(!is.na(cpg_id))

        colnames(donors_in_pool) <- c()
        print(donors_in_pool)

        new_file <- paste0(pools_dir, pool, ".tsv") # saves with the input pool name - NOT the tenk10k pool name
        if (save) {
            fwrite(donors_in_pool, new_file, sep = "\t") # create a file with the correct sequencing library name
        } else {
            print(donors_in_pool)
        }

        # get number of external ID's (donors) that are in each sequencing library
        n_donors_in_pool <- tenk10k_all_pool_info %>%
            filter(tenk10k_pool == {{ tenk10k_pool }}) %>%
            select(external_id) %>%
            nrow()

        return(data.frame(pool = pool, donors_in_pool = n_donors_in_pool))
    }

    # for each sample/sequencing library/pool:
    output <- all_samples_in_pool %>%
        pull(V1) %>% # 1. creates file containing the donors in each sequencig library
        purrr::map(\(x) get_donors_in_libraries(
            pool = x,
            pools_dir = pools_dir,
            save = save
        )) %>% # 2. counts how many donors there are and saves a table with total n donors for each of the sequencing libraries
        # purrr::map(
        #     \(x) get_n_donors_pool(pool = x)
        # ) %>%
        list_rbind() # make a df with donor numbers for all the samples contained in the cellranger outs

    if (save) {
        # print("saving...")
        output_path <- glue("{pool_list_path}/libraries_nsamples_{pool_list_suffix}.txt")
        write_tsv(output, output_path, col_names = FALSE)
    }
    return(output)
}


# PREVIOUS RUNS:
get_n_donors_all_pools(
    pool_list_suffix = "240501",
    pool_list_path = "data_processing",
    save = TRUE,
    pools_dir = pools_dir
)
get_n_donors_all_pools(
    pool_list_suffix = "240524",
    pool_list_path = "data_processing",
    save = TRUE,
    pools_dir = pools_dir
)
get_n_donors_all_pools(
    pool_list_suffix = "240223",
    pool_list_path = "data_processing",
    save = TRUE,
    pools_dir = pools_dir
)
# get_n_donors_all_pools(
#     pool_list_suffix = "231213",
#     pool_list_path = "data_processing",
#     save = TRUE,
#     pools_dir = pools_dir
# )
# get_n_donors_all_pools(
#     pool_list_suffix = "231214",
#     pool_list_path = "data_processing",
#     save = TRUE,
#     pools_dir = pools_dir
# )


# NOTE: ^ script does not account for the accidental (or intentional) maxi-pools S0034-37a and S0034-37b, need to manually create these? see /directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/Demuxafy/preprocessing/TOB_make_pool_sample_lists.R


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


################################################################
################### BIOHEART Accidental maxi pools #############
################################################################

# While analysing the original outputs of vireo (for demultiplexing), we observed something interesting about the pools from S0040 to S0047
# essentially, accidental maxi pools were generated by mistake during capture; the pools that were combined have to do with the physical location on the chip:
# S0040 + S0044, S0041 + S0045, S0042 + S0046, S0043 + S0047
# to make things more confusing, for S0040 and S0042 we actually only have the "_re" repeated sequencing, which is not affected by this

# find full explanation by Ellie at: https://docs.google.com/document/d/117F_3-dmLEsvgJr_6dRsZhg7lVZ8nSgduRLpajDBwCk/edit

# First, define the new expected list of donors for the maxi pools:
bh_file <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/TenK10K_BioHeart_pool_info_v5.csv"
bh_df <- read.csv(bh_file)

cpg_map_file <- "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/str_sample-sex-mapping_sample_karyotype_sex_mapping.csv"
cpg_map_df <- read.csv(cpg_map_file)

#### S0040 + S0044
tenk10k_pools <- c("S0040", "S0044")
sc_samples <- bh_df[bh_df$Pool %in% tenk10k_pools, "Sample.ID"]
cpg_samples <- cpg_map_df[cpg_map_df$external_id %in% sc_samples, "s"]
df <- data.frame(sample = cpg_samples)
new_file <- paste0(bioheart_samples_dir, "maxi_pool_S0044.tsv")
colnames(df) <- c()
fwrite(df, new_file, sep = "\t")

#### S0041 + S0045
tenk10k_pools <- c("S0041", "S0045")
sc_samples <- bh_df[bh_df$Pool %in% tenk10k_pools, "Sample.ID"]
cpg_samples <- cpg_map_df[cpg_map_df$external_id %in% sc_samples, "s"]
df <- data.frame(sample = cpg_samples)
new_file1 <- paste0(bioheart_samples_dir, "maxi_pool_S0041.tsv")
new_file2 <- paste0(bioheart_samples_dir, "maxi_pool_S0045.tsv")
colnames(df) <- c()
fwrite(df, new_file1, sep = "\t")
fwrite(df, new_file2, sep = "\t")

#### S0042 + S0046
tenk10k_pools <- c("S0042", "S0046")
sc_samples <- bh_df[bh_df$Pool %in% tenk10k_pools, "Sample.ID"]
cpg_samples <- cpg_map_df[cpg_map_df$external_id %in% sc_samples, "s"]
df <- data.frame(sample = cpg_samples)
new_file <- paste0(bioheart_samples_dir, "maxi_pool_S0046.tsv")
colnames(df) <- c()
fwrite(df, new_file, sep = "\t")

#### S0043 + S0047
tenk10k_pools <- c("S0043", "S0047")
sc_samples <- bh_df[bh_df$Pool %in% tenk10k_pools, "Sample.ID"]
cpg_samples <- cpg_map_df[cpg_map_df$external_id %in% sc_samples, "s"]
df <- data.frame(sample = cpg_samples)
new_file1 <- paste0(bioheart_samples_dir, "maxi_pool_S0043.tsv")
new_file2 <- paste0(bioheart_samples_dir, "maxi_pool_S0047.tsv")
colnames(df) <- c()
fwrite(df, new_file1, sep = "\t")
fwrite(df, new_file2, sep = "\t")

# Next, determine the numbers of expected  donors
output_file <- "/share/ScratchGeneral/anncuo/tenk10k/data_processing/libraries_nsamples_bioheart_accidental_maxipools.txt"

df <- data.frame(sample = paste0("S00", c(41, 43:47)))

# expected numbers (CPG IDs)
n <- nrow(fread(paste0(bioheart_samples_dir, "maxi_pool_S0044.tsv")))
df[df$sample %in% c("S0044"), "n_expected"] <- n
n <- nrow(fread(paste0(bioheart_samples_dir, "maxi_pool_S0041.tsv")))
df[df$sample %in% c("S0041", "S0045"), "n_expected"] <- n
n <- nrow(fread(paste0(bioheart_samples_dir, "maxi_pool_S0046.tsv")))
df[df$sample %in% c("S0046"), "n_expected"] <- n
n <- nrow(fread(paste0(bioheart_samples_dir, "maxi_pool_S0043.tsv")))
df[df$sample %in% c("S0043", "S0047"), "n_expected"] <- n

# all
n <- nrow(bh_df[bh_df$Pool %in% c("S0040", "S0044"), ])
df[df$sample %in% c("S0044"), "n_all"] <- n
n <- nrow(bh_df[bh_df$Pool %in% c("S0041", "S0045"), ])
df[df$sample %in% c("S0041", "S0045"), "n_all"] <- n
n <- nrow(bh_df[bh_df$Pool %in% c("S0042", "S0046"), ])
df[df$sample %in% c("S0046"), "n_all"] <- n
n <- nrow(bh_df[bh_df$Pool %in% c("S0043", "S0047"), ])
df[df$sample %in% c("S0043", "S0047"), "n_all"] <- n

# save
colnames(df) <- c()
fwrite(df, output_file, sep = "\t")


################################################################
##################### Old Multiome Pools #######################
################################################################

# Explanation from Ellie (https://centrepopgen.slack.com/archives/D033C82R1L5/p1706498219785809):

# Basically these four pools (Anna: S0056-59) were tested in 2021 for multiome under the name ES001_P1 through ES001_P4.
# We thought we'd be able to just thaw the existing pools and recapture, but on thawing there weren't enough cells of high enough quality.
# We captured what we could that day (which ended up being ES001_P1 (now called S0056) and a maxi-pool of ES001_P3 and ES001_P4 (now called S0058-59)
# knowing that the numbers would be low and I (=Ellie) would need to recreate these pools.
# When thawing the original PBMC vials again to repeat the pooling, some of the samples were insufficient to repeat, so didn't exist in the final iteration of the pool

# donor info recap (original vs new) at: https://docs.google.com/document/d/1pDB8argxXxOCsJ48qVAGZcfevNwAFiKFidDcFSkc6vs/edit

#### S0056

# what is called S0056 is an old pool orginally done for multiome, with 16 individuals
# it was then re-pooled (S0056_repeat) but there was only remaining material for 10 individuals (out of the 16)

S0056_donors <- c("CT_129", "CT_240", "CT_847", "CT_615", "CT_1069", "CT_1003", "CT_582", "CT_989", "CT_785", "CT_411")

S0056_original_donors <- c(
    "CT_985", "CT_129", "CT_127", "CT_240", "CT_847", "CT_615", "CT_1069", "CT_522", "CT_462",
    "CT_1003", "CT_685", "CT_582", "CT_85", "CT_989", "CT_785", "CT_411"
)

# note that CT_85 does not have a CPG ID correspective

sc_samples <- S0056_original_donors
cpg_samples <- cpg_map_df[cpg_map_df$external_id %in% sc_samples, "s"]
df <- data.frame(sample = cpg_samples)

new_file <- paste0(bioheart_samples_dir, "S0056.tsv")
colnames(df) <- c()
fwrite(df, new_file, sep = "\t")

#### next, the "S0058-59v2" maxipool

# S0058-59v2 is an old maxi pool of S0058 and S0059, orginally done for multiome
# S0058 "original" and S0059 "original" both had 16 donors
# S0058 and S0059 were then repeated, this time with only material for 7 and 12 donors respectively

S0058_donors <- c("CT_1036", "CT_137", "CT_164", "CT_76", "CT_170", "CT_294", "CT_788")

S0058_original_donors <- c(
    "CT_330", "CT_395", "CT_1036", "CT_633", "CT_78", "CT_137", "CT_164", "CT_722",
    "CT_76", "CT_170", "CT_714", "CT_112", "CT_315", "CT_294", "CT_788", "CT_1218"
)

S0059_donors <- c(
    "CT_207", "CT_116", "CT_421", "CT_1044", "CT_1073", "CT_470", "CT_398",
    "CT_824", "CT_156", "CT_134", "CT_107", "CT_341"
)

S0059_original_donors <- c(
    "CT_242", "CT_207", "CT_116", "CT_421", "CT_1044", "CT_1073", "CT_470", "CT_398",
    "CT_824", "CT_156", "CT_77", "CT_134", "CT_290", "CT_265", "CT_107", "CT_341"
)

sc_samples <- c(S0058_original_donors, S0059_original_donors)
cpg_samples <- cpg_map_df[cpg_map_df$external_id %in% sc_samples, "s"]
df <- data.frame(sample = cpg_samples)

new_file <- paste0(bioheart_samples_dir, "S0058-59v2.tsv")
colnames(df) <- c()
fwrite(df, new_file, sep = "\t")
