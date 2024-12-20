library(tidyverse)
library(glue)
library(data.table)

data_type <- "filtered_min1000genes_all_donors"

cell_metadata <- read_csv(glue("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/300_libraries/300_libraries_cell_metadata_{data_type}.csv")) %>%
    rename("barcode" = 1)

cell_metadata$individual %>%
    unique() %>%
    length()

indiv_list <- cell_metadata$individual %>%
    unique()

length(indiv_list)

cell_metadata$sequencing_library %>%
    unique() %>%
    length()


indiv_list_katie <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/temp/1924_qc_pass_samples.csv")

indiv_list_katie$s %>% length()

unmatched <- setdiff(
    indiv_list,
    indiv_list_katie$s
)

cell_metadata %>%
    filter(individual == unmatched) %>%
    group_by(individual, ct_id, sequencing_library) %>%
    count()



cell_metadata %>%
    colnames()

wgs_qc_df <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/temp/saige-qtl_tenk10k-genome-2-3-all_samples_annotated_by_filter.csv")


donors_missing <- c("CT_376", "CT_1037", "CT_1358", "CT_1568", "CT_1928", "CT_1936", "TOB1570", "TOB1730", "TOB1181", "TOB0901")

# CT_376 S0046	CPG247973 # abnormal cell type distribution
# CT_1037	S0052	CPG249177 # abnormal cell type distribution
# CT_1568	S0144	CPG251793 # abnormal cell type distribution
# CT_1358 # unknown - No record in https://docs.google.com/spreadsheets/d/1XcAYd3ssKTfdPquQ3QxtN8Mn18YObpiiuUecQFbz45Q/edit?gid=336263672#gid=336263672
# CT_1928 S0112	CPG254029 #  No cells demultiplexed by Vireo
# CT_1936 S0145	CPG254060 #  No cells demultiplexed by Vireo

# TOB1570 CPG309724 # abnormal cell type distribution
# TOB1730 CPG310938 # abnormal cell type distribution
# TOB1181 CPG312025 # abnormal cell type distribution
# TOB0901 CPG315986 # abnormal cell type distribution

search_ids <- c("CPG254029", "CPG254060")
search_ids <- c("CPG254029", "CPG254060")

search_ids %in% cell_metadata$individual

"CT_1358" %in% cell_metadata$ct_id
