library(data.table)

donor_dir = "/share/ScratchGeneral/anncuo/tenk10k/donor_info/"
bh_file = paste0(donor_dir,"TenK10K_BioHeart_pool_info.csv")
bh_df = read.csv(bh_file)

bioheart_samples_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/samples_in_pools_bioheart/"
tenk10k_pools = unique(bh_df$Pool)

for (tenk10k_pool in tenk10k_pools){
    df = data.frame(sample = bh_df[bh_df$Pool == tenk10k_pool,"Sample.ID"])
    new_file = paste0(bioheart_samples_dir, tenk10k_pool,".tsv")
    colnames(df) <- c()
    fwrite(df, new_file, sep="\t")
}
