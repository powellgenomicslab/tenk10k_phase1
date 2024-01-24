library(data.table)

# BioHEART donor information (which BioHEART IDs are in which pool)
donor_dir = "/share/ScratchGeneral/anncuo/tenk10k/donor_info/"
bh_file = paste0(donor_dir,"TenK10K_BioHeart_pool_info_v2.csv")
bh_df = read.csv(bh_file)

# Map between BioHEART IDs and CPG IDs (from Hope)
cpg_map_file = "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/data_processing/str_sample-sex-mapping_sample_karyotype_sex_mapping.csv"
cpg_map_df = read.csv(cpg_map_file)

# Output directory 
bioheart_samples_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/samples_in_pools_bioheart/"
tenk10k_pools = unique(bh_df$Pool)

# Loop over pools and create files
for (tenk10k_pool in tenk10k_pools){
    # CT IDs
    sc_samples = bh_df[bh_df$Pool == tenk10k_pool,"Sample.ID"]
    # corresponding CPG IDs
    cpg_samples = cpg_map_df[cpg_map_df$external_id %in% sc_samples, "s"]
    # build dataframe
    df = data.frame(sample = cpg_samples)
    # save
    new_file = paste0(bioheart_samples_dir, tenk10k_pool,".tsv")
    colnames(df) <- c()
    fwrite(df, new_file, sep="\t")
}

# for repeat samples just copy files over (S0040_repeat, S0042_repeat, S0103_re)

# for the one maxi pool (S0058-S0059v2), do this ad hoc:
tenk10k_pools = c("S0058","S0059")
sc_samples = bh_df[bh_df$Pool %in% tenk10k_pools,"Sample.ID"]
cpg_samples = cpg_map_df[cpg_map_df$external_id %in% sc_samples, "s"]
df = data.frame(sample = cpg_samples)
new_file = paste0(bioheart_samples_dir, "S0058-59v2.tsv")
colnames(df) <- c()
fwrite(df, new_file, sep="\t")
