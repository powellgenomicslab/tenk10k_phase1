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

###################

# save pool donor counts file
output_file = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/libraries_nsamples.txt"

cellranger_samples = list.files("/directflow/SCCGGroupShare/projects/data/experimental_data/projects/TenK10K/GencodeV44/","S")
simple_samples = cellranger_samples[!(cellranger_samples %in% cellranger_samples[grep("-",cellranger_samples)])]

# simple pools
for (sample in simple_samples){
    df_s = read.csv(paste0(samples_dir,sample,".tsv"),sep="\t", header = F)
    pool = gsub("a","",gsub("b","",gsub("c","",sample)))
    df[df$V1 == sample, "V2"] = nrow(df_s)
    df[df$V1 == sample, "V3"] = length(unique(pool_df[pool_df$Tenk10k_pool == pool, "TOB_ID"]$TOB_ID))
}

# add maxi pools
df[df$V1 == "S0021-24a", "V3"] = df[df$V1 == "S0021a", "V3"] + df[df$V1 == "S0022a", "V3"] + df[df$V1 == "S0023a", "V3"] + df[df$V1 == "S0024a", "V3"]
df[df$V1 == "S0021-24b", "V3"] = df[df$V1 == "S0021-24a", "V3"]
df[df$V1 == "S0025-28a", "V3"] = df[df$V1 == "S0025a", "V3"] + df[df$V1 == "S0026a", "V3"] + df[df$V1 == "S0027a", "V3"] + df[df$V1 == "S0028a", "V3"]
df[df$V1 == "S0025-28b", "V3"] = df[df$V1 == "S0025-28a", "V3"]
df[df$V1 == "S0030-33a", "V3"] = df[df$V1 == "S0030a", "V3"] + df[df$V1 == "S0031a", "V3"] + df[df$V1 == "S0032a", "V3"] + df[df$V1 == "S0033a", "V3"]
df[df$V1 == "S0030-33b", "V3"] = df[df$V1 == "S0030-33a", "V3"]

# save
colnames(df) <- c()
fwrite(df, output_file, sep="\t")
