library(tidyverse)
library(ggsci)

# This script can be sourced to get colour palettes, functions, and plot themes for consistent styling etc. for tenk10k studies


# ⬇️ Data getters ----

# use this to read in the latest metadata csv
# when the path changes, can just change this path so we don't need to update it in every script

get_latest_metadata <- function(
    # update these when the object path changes
    csv = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_cell_metadata_subset_filtered_reanalysed.csv") {
    metadata <- read_csv(csv) %>%
        rename("barcode" = 1) %>%
        mutate(cell_type = str_replace_all(wg2_scpred_prediction, "_", " "))
    return(metadata)
}

get_latest_umap <- function(
    # update these when the object path changes
    csv = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_harmony_umap_coords_filtered_reanalysed.csv") {
    metadata <- read_csv(csv) %>%
        rename("barcode" = 1)
    return(metadata)
}

# 🌈 colour palettes ----

# generate nested colors with material theme
red <- pal_material("red")(5)[2:5]
purple <- pal_material("purple")(4)[2:4]
pink <- pal_material("pink")(7)
deeppurple <- pal_material("deep-purple")(7)
indigo <- pal_material("indigo")(7)[]
blue <- pal_material("blue")(7)[]
lightblue <- pal_material("light-blue")(7)
teal <- pal_material("teal")(7)[]
lightgreen <- pal_material("light-green")(5)[2:5]
lime <- pal_material("lime")(7)[]
yellow <- pal_material("yellow")(7)[]
orange <- pal_material("orange")(7)[]
brown <- pal_material("brown")(4)[2:4]
grey <- pal_material("grey")(4)[2:4]

tenk_color_pal <- tribble(
    ~wg2_scpred_prediction, ~color,
    # Lymphoid
    ## B cells
    "B_intermediate", red[1],
    "B_memory", red[4],
    "B_naive", orange[4],
    "Plasmablast", orange[7],
    ## NK cells
    "NK", pink[5],
    "NK_CD56bright", purple[1],
    "NK_Proliferating", purple[3],
    # CD8 T cells
    "CD8_Naive", deeppurple[3],
    "CD8_Proliferating", deeppurple[6],
    "CD8_TCM", indigo[4],
    "CD8_TEM", indigo[7],
    # CD4 T cells
    "CD4_CTL", blue[2],
    "CD4_Naive", blue[4],
    "CD4_Proliferating", blue[6],
    "CD4_TCM", lightblue[1],
    "CD4_TEM", lightblue[4],
    "Treg", lightblue[6],
    #
    "dnT", grey[3],
    "gdT", grey[6],
    "ILC", teal[3],
    "MAIT", teal[5],
    # myeloid
    ## DC
    "pDC", lime[3],
    "cDC1", lime[7],
    "cDC2", lightgreen[3],
    "ASDC", lightgreen[4],
    ## Monocyte
    "CD14_Mono", yellow[2],
    "CD16_Mono", yellow[6],
    #
    "HSPC", brown[3],
    #
    # "Platelet", brown[1],
    # "Eryth", brown[2],
    #
    # "Doublet", grey[2],
) %>%
    mutate(cell_type = str_replace_all(wg2_scpred_prediction, "_", " "))

# ggplot themes ----

# define a ggplot theme that we can re-use across all the plots for consistent styling
# work in progress

theme_tenk10k <- function() {
    font <- "Helvetica"

    theme_classic() %+replace%
        theme(
            plot.title = element_blank(),
        )
}
