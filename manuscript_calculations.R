library(tidyverse)
library(data.table)

sc_metadata_all_donors <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_filtered_min1000genes_all_donors.csv")

## QC APPLIED TO THIS OBJECT:
# donor QC:
#       donors who withdrew consent are removed
#       The odd Celltype composition donors, and WGS-QC fails are all included
#       donors with any number of cells are included

# single-cell qc:
# single cell QC was performed with the following for consistency with Anna's analysis:
#       ["pct_counts_mt"] < 20
#       ["n_genes_by_counts"] > 1000
#       ["n_genes_by_counts"] < 10000
#       ["total_counts"] > 800


wgs_qc_fails <- read_csv(
    "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/saige-qtl_tenk10k-genome-2-3-eur_all_samples_to_drop.csv"
)

# from https://docs.google.com/spreadsheets/d/1NUqLIo0PSKbiFgNCCFpuAfyTRtwVciW7ceFreEt1Tsw/edit?gid=66312870#gid=66312870, filters = {diverse ancestry}
diverse_ancestry_only_indiv <- c(
    "CPG247437",
    "CPG247452",
    "CPG247460",
    "CPG247544",
    "CPG247684",
    "CPG247718",
    "CPG247775",
    "CPG247874",
    "CPG248104",
    "CPG248120",
    "CPG248203",
    "CPG248369",
    "CPG248401",
    "CPG248419",
    "CPG248435",
    "CPG248468",
    "CPG248591",
    "CPG248617",
    "CPG248724",
    "CPG248765",
    "CPG248856",
    "CPG248864",
    "CPG248906",
    "CPG249060",
    "CPG249136",
    "CPG249185",
    "CPG249227",
    "CPG249235",
    "CPG249243",
    "CPG249367",
    "CPG249409",
    "CPG249417",
    "CPG249482",
    "CPG249532",
    "CPG249599",
    "CPG249623",
    "CPG249664",
    "CPG249748",
    "CPG249839",
    "CPG249847",
    "CPG249987",
    "CPG250027",
    "CPG250050",
    "CPG250159",
    "CPG250183",
    "CPG250233",
    "CPG250266",
    "CPG250373",
    "CPG250480",
    "CPG250506",
    "CPG250597",
    "CPG250696",
    "CPG250860",
    "CPG251173",
    "CPG251199",
    "CPG251413",
    "CPG251454",
    "CPG251496",
    "CPG251520",
    "CPG251579",
    "CPG251629",
    "CPG251637",
    "CPG251835",
    "CPG251843",
    "CPG251926",
    "CPG251959",
    "CPG251991",
    "CPG252080",
    "CPG252114",
    "CPG252239",
    "CPG252387",
    "CPG252411",
    "CPG252437",
    "CPG252478",
    "CPG252676",
    "CPG252841",
    "CPG252858",
    "CPG252908",
    "CPG252932",
    "CPG252973",
    "CPG253146",
    "CPG253211",
    "CPG253245",
    "CPG253252",
    "CPG253435",
    "CPG253468",
    "CPG253567",
    "CPG253575",
    "CPG253591",
    "CPG253625",
    "CPG253823",
    "CPG253831",
    "CPG253922",
    "CPG253948",
    "CPG254326",
    "CPG254342",
    "CPG254466",
    "CPG254474",
    "CPG254573",
    "CPG254599",
    "CPG254607",
    "CPG254623",
    "CPG254672",
    "CPG254797",
    "CPG254870",
    "CPG254953",
    "CPG254979",
    "CPG255091",
    "CPG255299",
    "CPG255349",
    "CPG255380",
    "CPG255455",
    "CPG255588",
    "CPG255703",
    "CPG255810",
    "CPG256081",
    "CPG256115",
    "CPG256156",
    "CPG256164",
    "CPG256172",
    "CPG256230",
    "CPG256339",
    "CPG256420",
    "CPG305375",
    "CPG305557",
    "CPG305565",
    "CPG305649",
    "CPG305714",
    "CPG305748",
    "CPG305805",
    "CPG306092",
    "CPG306100",
    "CPG306134",
    "CPG306167",
    "CPG306175",
    "CPG313296",
    "CPG313676",
    "CPG315838",
    "CPG316059",
    "CPG316364",
    "CPG316513",
    "CPG317016",
    "CPG498782",
    "CPG498840",
    "CPG498873",
    "CPG498949",
    "CPG498980",
    "CPG499012",
    "CPG499129",
    "CPG499251",
    "CPG499384",
    "CPG499509",
    "CPG499533",
    "CPG499558",
    "CPG499574",
    "CPG499590",
    "CPG499624",
    "CPG499632",
    "CPG499731",
    "CPG499806",
    "CPG499905",
    "CPG499970",
    "CPG500041",
    "CPG500272",
    "CPG500280",
    "CPG500298",
    "CPG500322",
    "CPG500389",
    "CPG500629",
    "CPG500751",
    "CPG500769",
    "CPG500777",
    "CPG500801",
    "CPG500868",
    "CPG500884",
    "CPG500942",
    "CPG501007",
    "CPG501056",
    "CPG501106",
    "CPG501130",
    "CPG501148",
    "CPG501197",
    "CPG501213",
    "CPG501247",
    "CPG501387",
    "CPG501486",
    "CPG501528",
    "CPG501601",
    "CPG501627",
    "CPG501676",
    "CPG501742",
    "CPG501767",
    "CPG501874",
    "CPG501916",
    "CPG501940",
    "CPG507962",
    "CPG508051",
    "CPG508119",
    "CPG508176",
    "CPG508259"
 )

# get WGS qc faile, not including individuals excluded for diverse ancestry 
wgs_qc_fails_not_diverse_ancestry <- setdiff(wgs_qc_fails$s, diverse_ancestry_only_indiv)

sc_metadata_qc_filtered <- sc_metadata_all_donors %>%
    group_by(cpg_id) %>%
    filter(n() > 100) %>% # remove donors with less than 100 cells
    ungroup() %>%
    filter(!cpg_id %in% c(
        "CPG309724",
        "CPG310938",
        "CPG312025",
        "CPG315986",
        "CPG247973",
        "CPG249177",
        "CPG251793",
        "CPG252494",
        "CPG254169",
        "CPG254318",
        "CPG255760",
        "CPG249904"
    )) %>% # remove donors with abnormal cell type compositions
    filter(!cpg_id %in% wgs_qc_fails_not_diverse_ancestry) # remove wgs qc fails but keep those with non-european ancestry

# sanity checks
sc_metadata_qc_filtered %>%
    select(cpg_id) %>%
    n_distinct() # n donors matches the number we expect

# number of barcodes
sc_metadata_qc_filtered %>% nrow()
sc_metadata_qc_filtered %>%
    select(`...1`) %>%
    n_distinct()


# read in the age data
sample_meta <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/saige-qtl_tenk10k-genome-2-3-eur_input_files_241210_covariates_sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv")

sample_meta_cohorts_including_noneuropean <- sc_metadata_qc_filtered %>%
    select("cpg_id", "cohort") %>%
    distinct() %>% 
    left_join(sample_meta, by = c("cpg_id" = "sample_id"))

# # calculate cohort median age
medianage <- sample_meta_cohorts_including_noneuropean %>%
    group_by(cohort) %>%
    summarise(median_age = median(age))

sample_meta_cohorts_including_noneuropean
