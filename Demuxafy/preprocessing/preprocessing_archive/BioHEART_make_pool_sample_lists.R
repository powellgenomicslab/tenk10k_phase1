library(tidyverse)
library(data.table)
library(glue)

# BioHEART donor information (which BioHEART IDs are in which pool)
bioheart_samples_dir <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/samples_in_pools_bioheart/"

################################################################
################### Most BioHEART pools ########################
################################################################

# if this shows 0 donors in a pool check that the TenK10K_BioHeart_pool_info table has been updated, and check that the column names have not been changed
get_n_donors_all_pools <- function(pool_list_suffix, pool_list_path, save = FALSE, bioheart_samples_dir) {
  # get list of all samples in the sequencing library
  all_samples_in_pool <- fread(glue("{pool_list_path}/cellranger_outs_{pool_list_suffix}.txt"), header = FALSE)

  # creates a file for each sequencing library showing the donors in that sample pool
  get_donors_in_libraries <- function(sample, bioheart_samples_dir, save) {
    # read in the inputs
    # donor_dir <- "/share/ScratchGeneral/anncuo/tenk10k/donor_info/"
    # bh_file <- paste0(donor_dir, "TenK10K_BioHeart_pool_info_v4.csv")
    bh_file <- "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/TenK10K_BioHeart_pool_info_v5.csv"
    bh_df <- read.csv(bh_file)

    cpg_map_file <- "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/str_sample-sex-mapping_sample_karyotype_sex_mapping.csv"
    cpg_map_df <- read.csv(cpg_map_file)

    # strip off the annotations denoting re-sequenced samples
    tenk10k_pool <- str_remove(sample, pattern = "_re$|_repeat$|a$|b$|_re1$")
    print(paste0("Sample: ", sample, " Pool: ", tenk10k_pool))
    # use the tables to determine which donors (CPG ID's) are in each pool
    cpg_samples <- bh_df[bh_df$Pool == tenk10k_pool, "CPG_ID"]
    df <- data.frame(sample = cpg_samples)
    colnames(df) <- c()
    print(df)

    new_file <- paste0(bioheart_samples_dir, sample, ".tsv") # save to specified directory
    if (save) {
      fwrite(df, new_file, sep = "\t") # create a file with the correct sequencing library name
    }
  }

  # this function gets the number of donors for a given sample
  get_n_donors_sample <- function(sample) {
    # strip off the annotation denoting where the same sample has been re-sequenced
    # tenk10k_pool_name <- str_remove(sample, pattern = "_re$|_repeat$|a$|b$")
    donors_in_pool_file <- glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/demuxafy/samples_in_pools_bioheart/{sample}.tsv") # file generated in the above function

    # get number of rows (donors) that are in each sequencing library
    donors_in_pool <- fread(donors_in_pool_file, header = FALSE) %>%
      nrow()
    return(data.frame(sample = sample, donors_in_pool = donors_in_pool))
  }

  # for each sample/sequencing library/pool:
  output <- all_samples_in_pool %>%
    pull(V1) %>% # 1. creates file containing the donors in each sequencig library
    purrr::walk(\(x) get_donors_in_libraries(
      sample = x,
      bioheart_samples_dir = bioheart_samples_dir,
      save = save
    )) %>% # 2. counts how many donors there are and saves a table with total n donors for each of the sequencing libraries
    purrr::map(
      \(x) get_n_donors_sample(sample = x)
    ) %>%
    list_rbind() # make a df with donor numbers for all the samples contained in the cellranger outs

  if (save) {
    # print("saving...")
    output_path <- glue("{pool_list_path}/libraries_nsamples_{pool_list_suffix}.txt")
    write_tsv(output, output_path, col_names = FALSE)
  }
  return(output)
}

# get_n_donors_all_pools(
#   pool_list_suffix = "240214",
#   pool_list_path = "data_processing",
#   save = "F",
#   bioheart_samples_dir = bioheart_samples_dir
# )

# get_n_donors_all_pools(
#   pool_list_suffix = "240223",
#   pool_list_path = "data_processing",
#   save = TRUE,
#   bioheart_samples_dir = bioheart_samples_dir
# )

# get_n_donors_all_pools(
#   pool_list_suffix = "possible_sample_swaps",
#   pool_list_path = "data_processing",
#   save = TRUE,
#   bioheart_samples_dir = bioheart_samples_dir
# )

get_n_donors_all_pools(
  pool_list_suffix = "240501",
  pool_list_path = "data_processing",
  save = TRUE,
  bioheart_samples_dir = bioheart_samples_dir
)

get_n_donors_all_pools(
  pool_list_suffix = "240524",
  pool_list_path = "data_processing",
  save = TRUE,
  bioheart_samples_dir = bioheart_samples_dir
)


# for repeat samples just copy files over (S0040_repeat, S0042_repeat, S0083_re, S0088_re, S0094_re, S0103_re)


################################################################
################### Accidental maxi pools ######################
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
