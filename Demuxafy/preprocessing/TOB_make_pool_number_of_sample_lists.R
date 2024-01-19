# the files created here will be encoded as:
# library / sample | number of individuals from this library in VCF used | number of total individuals pooled and sequenced in this library

# load TOB TenK10K pool info
donor_info_folder = "/share/ScratchGeneral/anncuo/tenk10k/donor_info/"
pool_df = fread(paste0(donor_info_folder,"TOB_TenK10K_pool_info_v3.csv"))
colnames(pool_df) = gsub(" ","_",colnames(pool_df))

# load sample lists (created in TOB_make_pool_sample_lists.R)
samples_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/samples_in_pools_drews_names/"

#############################################################################
##########################     13 / 10 / 2023     ###########################
#############################################################################

Samples recap: all TOB, 64 libraries in total:

* Pools sequenced from S0016 to S0037 + 3 maxi pools (S0021-24, S0025-28, S0030-33)
* All single pools are sequenced three times, e.g., S0018a, S0018b, S0018c
  * missing here are S0016a, S0017a, S0019a, S0027a, S0037b and S0037c 
  * S0029 is missing completely (was skipped when pooling)
* All maxi pools are sequenced twice, e.g., S0021-24a, S0021-24b

############################

# libraries sequenced on Oct 13th, 2023 (hence the name 231013)
output_file = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/libraries_nsamples_231013.txt"

cellranger_samples = list.files("/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/","S")
df = data.frame(V1 = cellranger_samples)

# add expected and genotyped donor counts
for (sample in unique(df$V1)){
    df_s = read.csv(paste0(samples_dir, sample, ".tsv"), sep="\t", header = F)
    pool = gsub("a","",gsub("b","",gsub("c","",sample)))
    df[df$V1 == sample, "V2"] = nrow(df_s)
    df[df$V1 == sample, "V3"] = length(unique(pool_df[pool_df$Tenk10k_pool == pool, "TOB_ID"]$TOB_ID))
}

# add maxi pools (not directly linked to a single pool)
df[df$V1 == "S0021-24a", "V3"] = df[df$V1 == "S0021a", "V3"] + df[df$V1 == "S0022a", "V3"] + df[df$V1 == "S0023a", "V3"] + df[df$V1 == "S0024a", "V3"]
df[df$V1 == "S0021-24b", "V3"] = df[df$V1 == "S0021-24a", "V3"]
df[df$V1 == "S0025-28a", "V3"] = df[df$V1 == "S0025a", "V3"] + df[df$V1 == "S0026a", "V3"] + df[df$V1 == "S0027a", "V3"] + df[df$V1 == "S0028a", "V3"]
df[df$V1 == "S0025-28b", "V3"] = df[df$V1 == "S0025-28a", "V3"]
df[df$V1 == "S0030-33a", "V3"] = df[df$V1 == "S0030a", "V3"] + df[df$V1 == "S0031a", "V3"] + df[df$V1 == "S0032a", "V3"] + df[df$V1 == "S0033a", "V3"]
df[df$V1 == "S0030-33b", "V3"] = df[df$V1 == "S0030-33a", "V3"]

# save
colnames(df) <- c()
fwrite(df, output_file, sep="\t")


#############################################################################
##########################     13 / 12 / 2023     ###########################
#############################################################################

Samples recap: all TOB, 24 libraries in total:

* Pools sequenced from S0001 to S0008 
* All single pools are sequenced three times, e.g., S0008a, S0008b, S0008c
* No maxi pools are included

############################

# libraries sequenced on Dec 13th, 2023 (hence the name 231213)
output_file = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/libraries_nsamples_231213.txt"

# initialise data frame with sample names
cellranger_samples = list.files("/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231213_tenk10k_gencode44/cellranger_outs/","S")
df = data.frame(V1 = cellranger_samples)

# add expected and genotyped donor counts
for (sample in unique(df$V1)){
    df_s = read.csv(paste0(samples_dir, sample, ".tsv"), sep="\t", header = F)
    pool = gsub("a","",gsub("b","",gsub("c","",sample)))
    df[df$V1 == sample, "V2"] = nrow(df_s)
    df[df$V1 == sample, "V3"] = length(unique(pool_df[pool_df$Tenk10k_pool == pool, "TOB_ID"]$TOB_ID))
}

# save
colnames(df) <- c()
fwrite(df, output_file, sep="\t")

#############################################################################
##########################     14 / 12 / 2023     ###########################
#############################################################################

Samples recap: mix of TOB and BioHEART, 41 libraries in total:

TOB (24 libraries):

* Pools sequenced 
  * S0009 to S0014 
  * S0038 and S0039
  * 1 maxi pool (S0034-37)
* All single pools are sequenced three times, e.g., S0010a, S0010b, S0010c
  * except S0038 is sequenced four times (S0038a, S0038b, S0038c, S0038d) 
  * missing here are S0014b, S0014c, S0039c 
* The maxi pools is sequenced twice, i.e., S0034-37a, S0034-37b

BioHEART (17 libraries and pools):

* Pools sequenced
  * S0045 to S0053
  * S0060 to S0067
* All pools are sequenced just once

############################

# onto libraries sequenced on Dec 14th, 2023 (231214)
# this includes BioHEART individuals / pools, which will be handled differently.
output_file = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/libraries_nsamples_231214.txt"

cellranger_samples = list.files("/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/231214_tenk10k_gencode44/cellranger_outs/","S")
bioheart_pools = cellranger_samples[25:41]
df = data.frame(V1 = cellranger_samples[!(cellranger_samples %in% bioheart_pools)])

# add expected and genotyped donor counts
for (sample in unique(df$V1)){
    df_s = read.csv(paste0(samples_dir, sample, ".tsv"), sep="\t", header = F)
    pool = gsub("a","",gsub("b","",gsub("c","",gsub("d","",sample))))
    df[df$V1 == sample, "V2"] = nrow(df_s)
    df[df$V1 == sample, "V3"] = length(unique(pool_df[pool_df$Tenk10k_pool == pool, "TOB_ID"]$TOB_ID))
}

# before dealing with maxi pool, load older libraries
df0 = read.csv("/share/ScratchGeneral/anncuo/tenk10k/data_processing/libraries_nsamples_231013.txt",sep="\t", header=F)

df[df$V1 == "S0034-37a", "V3"] = df0[df0$V1 == "S0034a", "V3"] + df0[df0$V1 == "S0035a", "V3"] + df0[df0$V1 == "S0036a", "V3"] + df0[df0$V1 == "S0037a", "V3"]
df[df$V1 == "S0034-37b", "V3"] = df[df$V1 == "S0034-37a", "V3"]

# save
colnames(df) <- c()
fwrite(df, output_file, sep="\t")


#############################################################################
##########################     15 / 01 / 2024     ###########################
#############################################################################

Samples recap: mix of TOB and BioHEART, 25 libraries in total:

TOB (8 libraries):

* Pools sequenced 
  * entirety of S0015 (S0015a, S0015b, S0015c)
  * previously missing sequences of S0014 (S0014b, S0014c), S0016 (S0016a) and S0037 (S0037b, S0037c)

BioHEART (17 libraries and pools):

* Pools sequenced
  * S0040 to S0044
  * S0054 to S0056 (+ one maxi pool of S0058+S0059 together)
  * S0135
  * S0138 and S00139
  * S0141 to S0145
* All pools are sequenced just once


# onto libraries sequenced on Jan 15th, 2024 (240115)
# this includes BioHEART individuals / pools, which will be handled differently.
output_file = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/libraries_nsamples_240115.txt"

cellranger_samples = list.files("/directflow/GWCCGPipeline/projects/deliver/GIMR_GWCCG_230201_JOSPOW_10x_Tenk10k/240115_tenk10k_gencode44/cellranger_outs/","S")
bioheart_pools = cellranger_samples[9:25]
df = data.frame(V1 = cellranger_samples[!(cellranger_samples %in% bioheart_pools)])

# add expected and genotyped donor counts
for (sample in unique(df$V1)){
    df_s = read.csv(paste0(samples_dir, sample, ".tsv"), sep="\t", header = F)
    pool = gsub("a","",gsub("b","",gsub("c","",gsub("d","",sample))))
    df[df$V1 == sample, "V2"] = nrow(df_s)
    df[df$V1 == sample, "V3"] = length(unique(pool_df[pool_df$Tenk10k_pool == pool, "TOB_ID"]$TOB_ID))
}

# save
colnames(df) <- c()
fwrite(df, output_file, sep="\t")


