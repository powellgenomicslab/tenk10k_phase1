library(data.table)

# TenK10K pool info 
# matching TenK10K pools (S** form) to OneK1K pools (C** form) for all individuals (TOB-ID form)
tenk10k_pool_info = "/share/ScratchGeneral/anncuo/tenk10k/donor_info/TOB_TenK10K_pool_info_v3.csv"
tenk10k_pool_df = fread(tenk10k_pool_info)
colnames(tenk10k_pool_df) = gsub(" ","_",colnames(tenk10k_pool_df))

# OneK1K pool info
# matching OneK1K pool (C** form, here IMB_ID) to alternative pool / sample notation (1:77)
# also contains "GENOME1_ID" (?)
onek1k_pool_info = "/share/ScratchGeneral/anncuo/tenk10k/donor_info/OneK1K_TOB_pools.csv"
onek1k_pool_df = fread(onek1k_pool_info)

# Folder containing files for each OneK1K pool (sample notation, "OneK1K_scRNA_Sample{1:77}"
# each file contains the individuals (IID form) included in that pool 
drews_onek1k_donor_pool_info = "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/WG1-GoogleCloud-execution/Demultiplexing/data/individuals_per_pool/"

# Directory were saving samples in pools for TenK10K
# using Drew's files and notation
new_samples_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/samples_in_pools_drews_names/"

# loop over TenK10K pools
tenk10k_pools = unique(tenk10k_pool_df$Tenk10k_pool)
for (tenk10k_pool in tenk10k_pools){
    onek1k_pools = unique(tenk10k_pool_df[tenk10k_pool_df$Tenk10k_pool == tenk10k_pool,"Onek1k_pool"])$Onek1k_pool
    onek1k_samples = onek1k_pool_df[onek1k_pool_df$IMB_ID %in% onek1k_pools,"SAMPLE_NAME"]$SAMPLE_NAME
    drews_onek1k_donor_pool_file1 = paste0(drews_onek1k_donor_pool_info, onek1k_samples[1],".tsv")
    onek1k_donors_df1 = fread(drews_onek1k_donor_pool_file1)
    if (length(onek1k_samples)==1){
        df = onek1k_donors_df1
    }
    if (length(o_samples)==2){
        drews_onek1k_donor_pool_file2 = paste0(drews_onek1k_donor_pool_info, onek1k_samples[2],".tsv")
        onek1k_donors_df2 = fread(drews_onek1k_donor_pool_file2)
        df = rbind(onek1k_donors_df1,onek1k_donors_df2)
    }
    new_file = paste0(new_samples_dir, tenk10k_pool,".tsv")
    colnames(df) <- c()
    fwrite(df, new_file, sep="\t")
}
# then copy each file over for a,b,c replicates

# maxi pools

maxi_pool = "S0021-24a"
df1 = read.csv(paste0(new_samples_dir,"S0021.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(new_samples_dir,"S0022.tsv"),sep="\t", header = F)
df3 = read.csv(paste0(new_samples_dir,"S0023.tsv"),sep="\t", header = F)
df4 = read.csv(paste0(new_samples_dir,"S0024.tsv"),sep="\t", header = F)
df5 = rbind(df1,df2,df3,df4)
colnames(df5) <- c()
out_file = paste0(new_samples_dir, maxi_pool, ".tsv")
fwrite(df5, out_file, sep="\t")
# the copy for maxi_pool="S0021-24b"

maxi_pool = "S0025-28a"
df1 = read.csv(paste0(new_samples_dir,"S0025.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(new_samples_dir,"S0026.tsv"),sep="\t", header = F)
df3 = read.csv(paste0(new_samples_dir,"S0027.tsv"),sep="\t", header = F)
df4 = read.csv(paste0(new_samples_dir,"S0028.tsv"),sep="\t", header = F)
df5 = rbind(df1,df2,df3,df4)
colnames(df5) <- c()
out_file = paste0(new_samples_dir, maxi_pool, ".tsv")
fwrite(df5, out_file, sep="\t")
# the copy for maxi_pool="S0025-28b"

maxi_pool = "S0030-33a"
df1 = read.csv(paste0(new_samples_dir,"S0030.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(new_samples_dir,"S0031.tsv"),sep="\t", header = F)
df3 = read.csv(paste0(new_samples_dir,"S0032.tsv"),sep="\t", header = F)
df4 = read.csv(paste0(new_samples_dir,"S0033.tsv"),sep="\t", header = F)
df5 = rbind(df1,df2,df3,df4)
colnames(df5) <- c()
out_file = paste0(new_samples_dir, maxi_pool, ".tsv")
fwrite(df5, out_file, sep="\t")
# the copy for maxi_pool="S0030-33b"
