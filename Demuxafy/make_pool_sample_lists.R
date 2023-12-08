library(data.table)

# output directory
samples_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/demuxafy/samples_in_pools/"

# input directory
donor_info_folder = "/share/ScratchGeneral/anncuo/tenk10k/donor_info/"

# pool info
pool_df = fread(paste0(donor_info_folder,"TOB_TenK10K_pool_info.csv"))
colnames(pool_df) = gsub(" ","_",colnames(pool_df))

# donors recruited
donor_df = fread(paste0(donor_info_folder,"TOB_recruited_donors_info.csv"))

# samples we have genotypes for
df_samples_file = "/share/ScratchGeneral/anncuo/OneK1K/scrna-seq_grch38_association_files_OneK1K_CPG_IDs.tsv"
df_samples = read.csv(df_samples_file, sep="\t")

# loop over pools
for (pool in unique(pool_df$Tenk10k_pool)){
  tob_ids = pool_df[pool_df$Tenk10k_pool == pool, "TOB_ID"]$TOB_ID
  onek1k_ids = donor_df[donor_df$TOB_ID %in% tob_ids, "STUDY_ID"]$STUDY_ID
  onek1k_ids_in_geno = onek1k_ids[onek1k_ids %in% df_samples$OneK1K_ID]
  samples_df = data.frame(onek1k_ids_in_geno)
  colnames(samples_df) <- c() 
  out_file = paste0(samples_dir, pool, ".tsv")
  fwrite(samples_df, out_file, sep="\t")
}
# then copy each file over for a,b,c replicates

# maxi pools

maxi_pool = "S0021-24a"
df1 = read.csv(paste0(samples_dir,"S0021.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(samples_dir,"S0022.tsv"),sep="\t", header = F)
df1 = read.csv(paste0(samples_dir,"S0023.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(samples_dir,"S0024.tsv"),sep="\t", header = F)
df5 = rbind(df1,df2,df3,df4)
colnames(df5) <- c()
out_file = paste0(samples_dir, maxi_pool, ".tsv")
fwrite(df5, out_file, sep="\t")
# the copy for maxi_pool="S0021-24b"

maxi_pool = "S0025-28a"
df1 = read.csv(paste0(samples_dir,"S0025.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(samples_dir,"S0026.tsv"),sep="\t", header = F)
df1 = read.csv(paste0(samples_dir,"S0027.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(samples_dir,"S0028.tsv"),sep="\t", header = F)
df5 = rbind(df1,df2,df3,df4)
colnames(df5) <- c()
out_file = paste0(samples_dir, maxi_pool, ".tsv")
fwrite(df5, out_file, sep="\t")
# the copy for maxi_pool="S0025-28b"

maxi_pool = "S0030-33a"
df1 = read.csv(paste0(samples_dir,"S0030.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(samples_dir,"S0031.tsv"),sep="\t", header = F)
df1 = read.csv(paste0(samples_dir,"S0032.tsv"),sep="\t", header = F)
df2 = read.csv(paste0(samples_dir,"S0033.tsv"),sep="\t", header = F)
df5 = rbind(df1,df2,df3,df4)
colnames(df5) <- c()
out_file = paste0(samples_dir, maxi_pool, ".tsv")
fwrite(df5, out_file, sep="\t")
# the copy for maxi_pool="S0030-33b"
