# script by Drew Neavin to keep track of OneK1K donors

library(data.table)
library(tidyverse)
library(DBI)

datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/WG1-GoogleCloud-execution/Imputation/data/"
setwd(datadir)


### Read in Data ###
psam <- fread(paste0(datadir,"/plink/hg19_onek1k.psam"))
sc_IDs <- fread(paste0(datadir,"1k1k_ids_master.csv"), sep = ",")

con <- dbConnect(RSQLite::SQLite(), "/directflow/SCCGGroupShare/projects/data/experimental_data/PUBLISHING/OneK1K/Database/OneK1K_MetadataDB.db")

dbListTables(con)
#  [1] "CELL_METADATA"       "CELL_PCA"            "CELL_TYPE_IDS"      
#  [4] "CELL_UMAP"           "CHROMIUM_PREP_NOTES" "CHROMIUM_QC_METRICS"
#  [7] "DISEASES"            "DISEASES_DETAILED"   "DISEASE_TYPES"      
# [10] "DONOR_CELL_PREP"     "DONOR_IDS"           "DONOR_METADATA"     
# [13] "GENOTYPING_QC"       "MEDICATIONS"         "MEDICATIONS_SUMMARY"
# [16] "OTHER_MEDICATIONS"   "PATIENT_METADATA"    "POOL_IDS"           
# [19] "RECRUITED_DONORS"    "SMOKING_STATUS"     

data.table(dbGetQuery(con, "SELECT * FROM DONOR_CELL_PREP"))

data.table(dbGetQuery(con, "SELECT * FROM GENOTYPING_QC"))

data.table(dbGetQuery(con, "SELECT * FROM RECRUITED_DONORS"))

data.table(dbGetQuery(con, "SELECT * FROM DONOR_IDS"))

data.table(dbGetQuery(con, "SELECT * FROM PATIENT_METADATA"))

data.table(dbGetQuery(con, "SELECT * FROM CHROMIUM_PREP_NOTES"))

data.table(dbGetQuery(con, "SELECT * FROM CHROMIUM_QC_METRICS"))

data.table(dbGetQuery(con, "SELECT * FROM DONOR_METADATA"))

data.table(dbGetQuery(con, "SELECT * FROM POOL_IDS"))

data.table(dbGetQuery(con, "SELECT * FROM CELL_METADATA"))

data.table(dbGetQuery(con, "SELECT * FROM DISEASES"))

data.table(dbGetQuery(con, "SELECT * FROM DISEASES_DETAILED"))

data.table(dbGetQuery(con, "SELECT * FROM SMOKING_STATUS"))

data.table(dbGetQuery(con, "SELECT * FROM MEDICATIONS_SUMMARY"))




### IDentify repeated pools
donor_cell_prep <- data.table(dbGetQuery(con, "SELECT * FROM DONOR_CELL_PREP"))
pool_ids <- data.table(dbGetQuery(con, "SELECT * FROM POOL_IDS"))

unique(donor_cell_prep$PREP_NOTES)

pool_ids[IMB_ID == "C0199"]
donor_cell_prep[grep("C0199",PREP_NOTES)]
donor_cell_prep[IMB_ID == "C0199"]
pool_ids[IMB_ID == "C0235"]

donor_cell_prep[grep("C0199",PREP_NOTES)]$TOB_ID %in% donor_cell_prep[IMB_ID == "C0199"]$TOB_ID

pool_ids[IMB_ID == "C0225"]
donor_cell_prep[grep("C0225",PREP_NOTES)]
donor_cell_prep[IMB_ID == "C0225"]
pool_ids[IMB_ID == "C0236"]


### Look for twins
recruited_donors <- data.table(dbGetQuery(con, "SELECT * FROM RECRUITED_DONORS"))
genotyping_qc <- data.table(dbGetQuery(con, "SELECT * FROM GENOTYPING_QC"))
pool_ids <- data.table(dbGetQuery(con, "SELECT * FROM POOL_IDS"))
donor_cell_prep_named <- donor_cell_prep[pool_ids, on = "IMB_ID"]



unique(recruited_donors$Notes)
## Twins = TOB0792 & TOB0795


### Remove Twins
donor_cell_prep_named_subset <- donor_cell_prep_named[!(TOB_ID %in% c("TOB00792", "TOB-00795"))]


### Remove pools 40 and 66
table(donor_cell_prep_named_subset$SAMPLE_NAME)
donor_cell_prep_named_subset <- donor_cell_prep_named_subset[!(SAMPLE_NAME %in% c("OneK1K_scRNA_Sample40", "OneK1K_scRNA_Sample66"))]


### Remove immune disorder patient
donor_cell_prep_named_subset <- donor_cell_prep_named_subset[TOB_ID != "TOB-01641"]



### Intersect with recruited donors for IDs to match in genotyping_qc
donor_cell_prep_named_subset_IDs <- unique(recruited_donors[,c("STUDY_ID", "FID", "IID", "TOB_ID")])[donor_cell_prep_named_subset, on = "TOB_ID"]


genotyping_qc_ID_matched <- donor_cell_prep_named_subset_IDs[,c("FID", "IID", "TOB_ID")][genotyping_qc, on = c("IID", "FID")]
genotyping_qc_ID_matched <- unique(genotyping_qc_ID_matched)


genotyping_qc_ID_matched <- genotyping_qc_ID_matched[!(is.na(TOB_ID))]

genotyping_qc_ID_matched[STUDY_ID == "5_5"]


### Check for duplicates
for (iid in genotyping_qc_ID_matched[duplicated(TOB_ID)]$IID){
    print(genotyping_qc_ID_matched[grep(iid, IID)])
}


### TOB-01051 was mislabeled and seems to be in the list twice as a result -> remove the row with this label from the recruited_donors1023
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[TOB_ID != "TOB-01051"]

### TOB1638 has two genotyped samples: 155_155 and 1086_155. Previous demultiplexing indicates 155_155 was used and had lots of cells, remove the other duplicate genotyping, also genotyping quality metrics are better
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "1086_155"]


### TOB-01639 has two genotyped samples: 156_156 and 1087_156. Previous demultiplexing indicates 156_156 was used and had lots of cells, remove the other duplicate genotyping, also genotyping quality metrics are better
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "1087_156"]


### TOB-01660 has two genotyped samples: 157_157 and 1088_157. Previous demultiplexing indicates doesn't have either (wasn't provided as an individual in the pool, possibly failed genotyping QC/imputation?), genotyping quality metrics are better for 1088_157
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "157_157"]


### TOB-01757 has two genotyped samples: 253_254 and 1090_254. Previous demultiplexing indicates doesn't have either (wasn't provided as an individual in the pool, possibly failed genotyping QC/imputation?), genotyping quality metrics are better for 1090_254
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "253_254"]


### TOB1759 has two genotyped samples: 255_256 and 1092_256. Previous demultiplexing indicates 255_256 was used and had lots of cells, remove the other duplicate genotyping, also genotyping quality metrics are better
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "1092_256"]


### TOB1361 has two genotyped samples: 505_506 and 1093_506. Previous demultiplexing indicates 505_506 was used and had lots of cells, but genotyping metrics are better for the other one, so will use 1093_506
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "505_506"]


### TOB1410 has two genotyped samples: 553_554 and 1094_554. Previous demultiplexing indicates 505_506 was used and had lots of cells, but genotyping metrics are better for the other one, so will use 1094_554
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "553_554"]


### TOB0967 has two genotyped samples: 844_845 and 1096_845. Previous demultiplexing doesn't have either but 844_845 has worse metrics so will use 1096_845
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "844_845"]


### TOB0968 has two genotyped samples: 845_846 and 1097_846. Previous demultiplexing used 845_846 and had lots of cells and better qc metrics
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "1097_846"]


### TOB0969 has two genotyped samples: 846_847 and 1098_847. Previous demultiplexing used 846_847 and had lots of cells and better qc metrics
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "1098_847"]


### TOB0970 has two genotyped samples: 847_848 and 1099_848. Previous demultiplexing used 847_848 and had lots of cells and better qc metrics
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "1099_848"]


### TOB1081 has two genotyped samples: 913_914 and 1100_914. Previous demultiplexing has neither but 913_914 has bad metrics, so will use 1100_914
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "913_914"]


### TOB1084 has two genotyped samples: 914_915 and 1101_915. Previous demultiplexing had 914_915 but has bad metrics, so will use 1101_915
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "914_915"]


### TOB1806 has two genotyped samples: 931_932 and 1102_932. Previous demultiplexing doesn't have either 931_932 but has bad metrics, so will use 1102_932
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "931_932"]


### TOB1800 has two genotyped samples: 925_926 and 1103_926. Previous demultiplexing doesn't have either 925_926 but has bad metrics, so will use 1103_926
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "925_926"]


### TOB1089 has two genotyped samples: 919_920 and 1104_920. Previous demultiplexing doesn't have either 919_920 but has bad metrics, so will use 1104_920
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "919_920"]

### TOB1745 has two genotyped samples: 241_242 and 1089_242.1089_242 is better quality so will use that one
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "241_242"]

### TOB1745 has two genotyped samples: 565_566 and 1095_566. 1095_566 is better quality so will use that one
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "565_566"]

### TOB-01049 has two genotyped samples: 1022_1023 and 1023_1024. 1022_1023 has better data quality so use that one
genotyping_qc_ID_matched <- genotyping_qc_ID_matched[STUDY_ID != "1023_1024"]




### Check for duplicates
donor_cell_prep_named_subset[TOB_ID %in% donor_cell_prep_named_subset[duplicated(TOB_ID)]$TOB_ID]


## Two duplicates: TOB-01002 and TOB-01003 are in pools 19 and 38




recruited_donors1023 <- recruited_donors[PatientID %in% sc_IDs$PatientID]

### For some reason the individuals that should be in 76 and 77 (same as 40 and 66) are not all listed in the file (just one for each)
## Remove the 76 and 77 lone donors
# recruited_donors1023 <- recruited_donors1023[!(POOL_ID %in% c(76, 77))]

## Update pool IDs for 40 and 66 to 76 and 77
recruited_donors1023$POOL_ID <- gsub(40, 76, recruited_donors1023$POOL_ID) %>%
                                    gsub(66, 77, .)


### Remove Twins
recruited_donors1023 <- recruited_donors1023[!(PatientID %in% c("TOB0792", "TOB0795"))]


recruited_donors[TOB_ID %in% donor_cell_prep[grep("C0199",PREP_NOTES)]$TOB_ID]
recruited_donors[TOB_ID %in% donor_cell_prep[IMB_ID == "C0235"]$TOB_ID]

recruited_donors[TOB_ID %in% donor_cell_prep[grep("C0225",PREP_NOTES)]$TOB_ID]
recruited_donors[TOB_ID %in% donor_cell_prep[IMB_ID == "C0236"]$TOB_ID]




### Check for duplicates
genotyping_qc_1023 <- genotyping_qc
recruited_donors1023[PatientID %in% recruited_donors1023[duplicated(PatientID)]$PatientID]

for (iid in recruited_donors1023[duplicated(PatientID)]$IID){
    print(genotyping_qc_1023[grep(iid, IID)])
    print(recruited_donors1023[IID == iid])
}

### TOB-01051 was mislabeled and seems to be in the list twice as a result -> remove the row with this label from the recruited_donors1023
recruited_donors1023 <- recruited_donors1023[TOB_ID != "TOB-01051"]


### 869_870 everything is the same except DONOR_INDEX and only one genotyping_qc_1023 so remove this row
recruited_donors1023 <- recruited_donors1023[DONOR_INDEX != "873"]


### 870_871 has 4 rows between two pools. Weird note: If female, then TOB-01856, get the genotype id 967. But when check previous demultiplexing, can see this individual is truly in those pools but two rows per pool -> delete three of them 
recruited_donors1023 <- recruited_donors1023[DONOR_INDEX != "876"]
recruited_donors1023 <- recruited_donors1023[DONOR_INDEX != "877"]
recruited_donors1023 <- recruited_donors1023[DONOR_INDEX != "875"]


### TOB1638 has two genotyped samples: 155_155 and 1086_155. Previous demultiplexing indicates 155_155 was used and had lots of cells, remove the other duplicate genotyping, also genotyping quality metrics are better
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "1086_155"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "1086_155"]


### TOB-01639 has two genotyped samples: 156_156 and 1087_156. Previous demultiplexing indicates 156_156 was used and had lots of cells, remove the other duplicate genotyping, also genotyping quality metrics are better
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "1087_156"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "1087_156"]


### TOB-01660 has two genotyped samples: 157_157 and 1088_157. Previous demultiplexing indicates doesn't have either (wasn't provided as an individual in the pool, possibly failed genotyping QC/imputation?), genotyping quality metrics are better for 1088_157
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "157_157"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "157_157"]


### TOB-01757 has two genotyped samples: 253_254 and 1090_254. Previous demultiplexing indicates doesn't have either (wasn't provided as an individual in the pool, possibly failed genotyping QC/imputation?), genotyping quality metrics are better for 1090_254
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "253_254"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "253_254"]


### TOB1759 has two genotyped samples: 255_256 and 1092_256. Previous demultiplexing indicates 255_256 was used and had lots of cells, remove the other duplicate genotyping, also genotyping quality metrics are better
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "1092_256"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "1092_256"]


### TOB1361 has two genotyped samples: 505_506 and 1093_506. Previous demultiplexing indicates 505_506 was used and had lots of cells, but genotyping metrics are better for the other one, so will use 1093_506
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "505_506"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "505_506"]


### TOB1410 has two genotyped samples: 553_554 and 1094_554. Previous demultiplexing indicates 505_506 was used and had lots of cells, but genotyping metrics are better for the other one, so will use 1094_554
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "553_554"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "553_554"]


### TOB0967 has two genotyped samples: 844_845 and 1096_845. Previous demultiplexing doesn't have either but 844_845 has worse metrics so will use 1096_845
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "844_845"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "844_845"]


### TOB0968 has two genotyped samples: 845_846 and 1097_846. Previous demultiplexing used 845_846 and had lots of cells and better qc metrics
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "1097_846"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "1097_846"]


### TOB0969 has two genotyped samples: 846_847 and 1098_847. Previous demultiplexing used 846_847 and had lots of cells and better qc metrics
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "1098_847"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "1098_847"]


### TOB0970 has two genotyped samples: 847_848 and 1099_848. Previous demultiplexing used 847_848 and had lots of cells and better qc metrics
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "1099_848"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "1099_848"]


### TOB1081 has two genotyped samples: 913_914 and 1100_914. Previous demultiplexing has neither but 913_914 has bad metrics, so will use 1100_914
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "913_914"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "913_914"]


### TOB1084 has two genotyped samples: 914_915 and 1101_915. Previous demultiplexing had 914_915 but has bad metrics, so will use 1101_915
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "914_915"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "914_915"]


### TOB1806 has two genotyped samples: 931_932 and 1102_932. Previous demultiplexing doesn't have either 931_932 but has bad metrics, so will use 1102_932
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "931_932"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "931_932"]


### TOB1800 has two genotyped samples: 925_926 and 1103_926. Previous demultiplexing doesn't have either 925_926 but has bad metrics, so will use 1103_926
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "925_926"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "925_926"]


### TOB1089 has two genotyped samples: 919_920 and 1104_920. Previous demultiplexing doesn't have either 919_920 but has bad metrics, so will use 1104_920
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "919_920"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "919_920"]

### TOB1745 has two genotyped samples: 241_242 and 1089_242.1089_242 is better quality so will use that one
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "241_242"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "241_242"]

### TOB1745 has two genotyped samples: 565_566 and 1095_566. 1095_566 is better quality so will use that one 
# note by Anna: I think this is TOB1421 instead (not TOB1745), and I do not see duplicates
recruited_donors1023 <- recruited_donors1023[STUDY_ID != "565_566"]
genotyping_qc_1023 <- genotyping_qc_1023[STUDY_ID != "565_566"]


### TOB1696 has two gentries for pool 77
recruited_donors1023 <- recruited_donors1023[DONOR_INDEX != "192"]


### TOB1420 has two gentries for pool 76
recruited_donors1023 <- recruited_donors1023[DONOR_INDEX != "566"]


## TOB-01049 has two gentries for pool 16
recruited_donors1023 <- recruited_donors1023[DONOR_INDEX != "1030"]


genotyping_qc_1023


genotyping_qc_ID_matched
recruited_donors1023


any(duplicated(genotyping_qc_ID_matched$TOB_ID))
any(duplicated(genotyping_qc_ID_matched$IID))
any(duplicated(genotyping_qc_ID_matched$STUDY_ID))


any(duplicated(recruited_donors1023$TOB_ID))
any(duplicated(recruited_donors1023$IID))
any(duplicated(recruited_donors1023$STUDY_ID))

psam_subset <- psam[recruited_donors1023[,c("IID","FID", "PatientID")], on = c("IID", "#FID" = "FID")]


genotyping_qc_ID_matched_w_pool <- donor_cell_prep_named[genotyping_qc_ID_matched, on = c("TOB_ID")]
genotyping_qc_ID_matched_w_pool <- genotyping_qc_ID_matched_w_pool[!(SAMPLE_NAME %in% c("OneK1K_scRNA_Sample40", "OneK1K_scRNA_Sample66"))]

all(genotyping_qc_ID_matched_w_pool$STUDY_ID %in% recruited_donors1023$STUDY_ID)
genotyping_qc_ID_matched_w_pool[!(genotyping_qc_ID_matched_w_pool$STUDY_ID %in% recruited_donors1023$STUDY_ID)]
donor_cell_prep_named[TOB_ID %in% genotyping_qc_ID_matched_w_pool[!(genotyping_qc_ID_matched_w_pool$STUDY_ID %in% recruited_donors1023$STUDY_ID)]$TOB_ID]



all(recruited_donors1023$STUDY_ID %in% genotyping_qc_ID_matched$STUDY_ID)
recruited_donors1023[!(recruited_donors1023$STUDY_ID %in% genotyping_qc_ID_matched$STUDY_ID)]


genotyping_qc_ID_matched_w_pool[SAMPLE_NAME == "OneK1K_scRNA_Sample30"]
recruited_donors1023[POOL_ID == "30"]
genotyping_qc_ID_matched_w_pool[TOB_ID == "TOB-01239"]




recruited_donors1023[TOB_ID == recruited_donors1023[duplicated(TOB_ID)]$TOB_ID]
genotyping_qc_ID_matched[TOB_ID == "TOB-01049"]



unique(recruited_donors1023$Notes)

recruited_donors1023[Notes == "Sibling of TOB-00958"] ### These are siblings but in different pools so I will leave it so that WG3 can select what they want to do
recruited_donors1023[Notes == "Sibling of TOB-00801"] ### These are siblings but in different pools so I will leave it so that WG3 can select what they want to do



##### Update psam because it appears that some individuals have the same IID but different FID and so have _2
psam$IID <- gsub("_[0-9]$", "", psam$IID)

recruited_donors1023$IID <- as.character(recruited_donors1023$IID)

psam_subset <- psam[recruited_donors1023[,c("IID","FID", "PatientID")], on = c("IID", "#FID" = "FID")]


##### check if unique patient IDs
any(duplicated(psam_subset$PatientID))
any(duplicated(psam_subset$IID))




### Subset the clin_data for needed columns and updated to correct nomenclature for consortium ###
patient_metadata <- data.table(dbGetQuery(con, "SELECT * FROM PATIENT_METADATA"))
patient_metadata[, c("#FID", "IID") := tstrsplit(STUDY_ID, "_", fixed=TRUE)]

psam_subset$`#FID` <- as.character(psam_subset$`#FID`)



### Check if sex provided in the two tables match
psam_subset_updated <- patient_metadata[,c("PatientID", "AGE", "SEX")][psam_subset, on = c("PatientID")]

which(ifelse(psam_subset_updated$SEX == "Female" & psam_subset_updated$i.SEX == 1 | psam_subset_updated$SEX == "Male" & psam_subset_updated$i.SEX == 2, "Mismatch", "Match") == "Mismatch") 
## All match or are NA

### Remake table
psam_subset_updated <- patient_metadata[,c("PatientID", "AGE")][psam_subset, on = c("PatientID")]
colnames(psam_subset_updated) <- gsub("AGE", "age", colnames(psam_subset_updated))

psam_subset_updated$Provided_Ancestry <- "EUR"
psam_subset_updated$genotyping_platform <- "IlluminaInfiniumGlobalScreeningArray"
psam_subset_updated$array_available <- "Y"
psam_subset_updated$wgs_available <- "N"
psam_subset_updated$wes_available <- "N"
psam_subset_updated$age_range <- floor(psam_subset_updated$age/10)*10
psam_subset_updated$Study <- "OneK1K"
psam_subset_updated$smoking_status <- "NONE"
psam_subset_updated$hormonal_contraception_use_currently <- "NONE"
psam_subset_updated$menopause <- "NONE"
psam_subset_updated$pregnancy_status <- "NONE"



psam_final <- psam_subset_updated[,c("#FID","IID","PAT","MAT","SEX","Provided_Ancestry","genotyping_platform","array_available","wgs_available","wes_available","age","age_range","Study","smoking_status","hormonal_contraception_use_currently","menopause","pregnancy_status")]



### *** Write the list of individuals that will be needed to be subset from the data 
fwrite(psam_final[,c("#FID", "IID")], paste0(datadir, "samples2keep.txt"), sep = "\t")

### After using plink2 to subset in make_ped_files.sh, read in to check
psam_subset_plink <- fread(paste0(datadir, "plink_subset/hg19_onek1k_old.psam"))

all(psam_subset_plink$`#FID` == psam_final$`#FID`)
all(psam_subset_plink$`IID` == psam_final$`IID`)

### TRUE !!!


# psam_final[psam_final == "NaN"] = "NONE"
# psam_final[age == "NONE"]$age <- "NA"
# psam_final[age_range == "NONE"]$age_range <- "NA"

psam_final[73,]

fwrite(psam_final, paste0(datadir,"plink_subset/hg19_onek1k.psam"), sep = "\t", na = "NA", quote = FALSE)




##### Write files for demultiplexing and classification #####
### file with pools so know which to upload
demultiplexing_datadir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/WG1-GoogleCloud-execution/Demultiplexing/data/"
dir.create(demultiplexing_datadir)

fwrite(data.table(Pool =sort(paste0("OneK1K_scRNA_Sample",unique(recruited_donors1023[POOL_ID != "NA"]$POOL_ID)))), paste0(demultiplexing_datadir, "/pool_list.tsv"), col.names = FALSE)


### individuals per pool files
indiv_pool_dir <- paste0(demultiplexing_datadir, "individuals_per_pool/")
dir.create(indiv_pool_dir)

lapply(sort(unique(recruited_donors1023[POOL_ID != "NA"]$POOL_ID)), function(pool){
    print(paste0("OneK1K_scRNA_Sample", pool))
    print(length(recruited_donors1023[POOL_ID == pool]$STUDY_ID))
    # print(data.table(Individual = recruited_donors1023[POOL_ID == pool]$STUDY_ID))
    fwrite(data.table(Individual = recruited_donors1023[POOL_ID == pool]$STUDY_ID), paste0(indiv_pool_dir, "OneK1K_scRNA_Sample", pool, ".tsv"), col.names = FALSE)
})


### Sample dataframe that has the N indiviudlas per pool
pool_ids <- data.table(dbGetQuery(con, "SELECT * FROM POOL_IDS"))

donor_cell_prep_named <- donor_cell_prep[pool_ids, on = "IMB_ID"]



donor_cell_prep_named_table <- data.table(table(donor_cell_prep_named$SAMPLE_NAME))
colnames(donor_cell_prep_named_table) <- c("V1", "N_prep")

recruited_donors_table <- data.table(table(recruited_donors$POOL_ID))
recruited_donors_table$V1 <- paste0("OneK1K_scRNA_Sample", recruited_donors_table$V1)
colnames(recruited_donors_table) <- c("V1", "N_recruited")

sc_IDs_table <- data.table(table(sc_IDs$ChromiumSampleName))
colnames(recruited_donors_table) <- c("V1", "N_sc")


combined <- donor_cell_prep_named_table[recruited_donors_table, on = "V1"]
combined2 <- combined[sc_IDs_table, on = "V1"]
combined$difference <- as.numeric(combined$N_prep) - as.numeric(combined$N_recruited)


donor_cell_prep_named[SAMPLE_NAME == "OneK1K_scRNA_Sample62"]
recruited_donors[POOL_ID == 62]


donor_cell_prep_named[TOB_ID == "TOB-01590"]
recruited_donors[TOB_ID == "TOB-01590"]


recruited_donors[POOL_ID == "NA"]
donor_cell_prep_named[TOB_ID == "TOB-01555"]


donor_cell_prep_named[SAMPLE_NAME == "OneK1K_scRNA_Sample70"]
recruited_donors[POOL_ID == 70]
