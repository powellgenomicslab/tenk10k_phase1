library(tidyverse)
library(ggsci)

# This script can be sourced to get colour palettes, functions, and plot themes for consistent styling etc. for tenk10k studies

# ⬇️ Data getters ----

# use this to read in the latest metadata csv
# when the path changes, can just change this path so we don't need to update it in every script

get_latest_metadata <- function(
    # update these when the object path changes
    csv = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/240_libraries/240_libraries_cell_metadata_subset_filtered_reanalysed.csv") {
    # add in the clean cell type names for plotting
    metadata <- read_csv(csv) %>%
        rename("barcode" = 1) %>%
        mutate(
            cell_type = str_replace_all(wg2_scpred_prediction, "_", " "),
            major_cell_type = case_when(
                wg2_scpred_prediction %in% c(
                    "B_intermediate",
                    "B_memory",
                    "B_naive",
                    "Plasmablast"
                ) ~ "B",
                wg2_scpred_prediction %in% c(
                    "NK",
                    "NK_CD56bright",
                    "NK_Proliferating"
                ) ~ "NK",
                wg2_scpred_prediction %in% c(
                    "CD8_Naive",
                    "CD8_Proliferating",
                    "CD8_TCM",
                    "CD8_TEM"
                ) ~ "CD8 T",
                wg2_scpred_prediction %in% c(
                    "CD4_CTL",
                    "CD4_Naive",
                    "CD4_Proliferating",
                    "CD4_TCM",
                    "CD4_TEM",
                    "Treg"
                ) ~ "CD4 T",
                wg2_scpred_prediction %in% c(
                    "dnT",
                    "gdT",
                    "ILC",
                    "MAIT"
                ) ~ "Unconventional T",
                wg2_scpred_prediction %in% c(
                    "pDC",
                    "cDC1",
                    "cDC2",
                    "ASDC"
                ) ~ "Dendritic",
                wg2_scpred_prediction %in% c(
                    "CD14_Mono",
                    "CD16_Mono"
                ) ~ "Monocyte",
                wg2_scpred_prediction %in% c(
                    "HSPC"
                ) ~ "Other"
            ) %>% fct_relevel(c("B", "NK", "CD8 T", "CD4 T", "Unconventional T", "Dendritic", "Monocyte", "Other"))
        )
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
red <- pal_material("red")(10)
purple <- pal_material("purple")(10)
pink <- pal_material("pink")(10)
deeppurple <- pal_material("deep-purple")(10)
indigo <- pal_material("indigo")(10)
blue <- pal_material("blue")(10)
lightblue <- pal_material("light-blue")(10)
teal <- pal_material("teal")(10)
green <- pal_material("green")(10)
lime <- pal_material("lime")(10)
yellow <- pal_material("yellow")(10)
orange <- pal_material("orange")(10)
brown <- pal_material("brown")(10)
grey <- pal_material("grey")(10)

# NOTE: the "color" column is used in downstream plotting, modifying this will update the colour scheme for all downstream plots

tenk_color_pal <- tribble(
    ~wg2_scpred_prediction, ~color_material, ~color_biorender1, ~color,
    # Lymphoid
    # CD4 T cells
    "CD4_TCM", blue[10], "#AB728A", "#9C85C3",
    "CD4_Naive", blue[8], "#83576A", "#605278",
    "CD4_TEM", blue[6], "#C4839E", "#8773A8",
    "CD4_CTL", blue[4], "#DA91B0", "#B197DD",
    "Treg", blue[3], "#F4BFE0", "#DFB4EC",
    "CD4_Proliferating", blue[2], "#DA91B0", "#D4ACE1",
    # CD8 T cells
    "CD8_TEM", deeppurple[10], "#749C93", "#6A81B5",
    "CD8_Naive", deeppurple[7], "#57736E", "#5D6E8D",
    "CD8_TCM", deeppurple[5], "#70B1A4", "#A7B4D3",
    "CD8_Proliferating", deeppurple[2], "#8BD1BF", "#CDDAE7",
    ## NK cells
    "NK", pink[9], "#708A5F", "#708A5F",
    "NK_CD56bright", pink[6], "#96B980", "#96B981",
    "NK_Proliferating", pink[3], "#BAE49E", "#BAE49E",
    ## B cells
    "B_naive", yellow[10], "#EEE0B6", "#EEE0B7",
    "B_intermediate", yellow[8], "#F8ED7D", "#F8ED7D",
    "B_memory", yellow[6], "#E2D872", "#E2D872",
    "Plasmablast", yellow[4], "#B5AE5C", "#B5AE5C",
    # Unconventional T
    "gdT", lime[10], "#A7B4D3", "#DE94B4",
    "MAIT", lime[8], "#5D6E8D", "#B0758E",
    "dnT", lime[6], "#6A81B5", "#EE9FC1",
    "ILC", lime[4], "#CDDAE7", "#BB7D97",

    ## Monocyte
    # myeloid
    "CD14_Mono", green[8], "#F3AF84", "#F3AF84",
    "CD16_Mono", green[4], "#DD9F78", "#DD9F78",
    ## DC
    "cDC2", brown[8], "#9C4C49", "#9C4C49",
    "pDC", brown[6], "#ED7470", "#ED7470",
    "cDC1", brown[4], "#6E3634", "#6E3634",
    "ASDC", brown[2], "#CF6662", "#CF6662",
    # Other
    "HSPC", grey[5], "#A9A9A9", "#A9A9A9",
    #
    # "Platelet", brown[1],
    # "Eryth", brown[2],
    #
    # "Doublet", grey[2],
) %>%
    mutate(
        cell_type = str_replace_all(wg2_scpred_prediction, "_", " "),
        major_cell_type = case_when(
            wg2_scpred_prediction %in% c(
                "B_intermediate",
                "B_memory",
                "B_naive",
                "Plasmablast"
            ) ~ "B",
            wg2_scpred_prediction %in% c(
                "NK",
                "NK_CD56bright",
                "NK_Proliferating"
            ) ~ "NK",
            wg2_scpred_prediction %in% c(
                "CD8_Naive",
                "CD8_Proliferating",
                "CD8_TCM",
                "CD8_TEM"
            ) ~ "CD8 T",
            wg2_scpred_prediction %in% c(
                "CD4_CTL",
                "CD4_Naive",
                "CD4_Proliferating",
                "CD4_TCM",
                "CD4_TEM",
                "Treg"
            ) ~ "CD4 T",
            wg2_scpred_prediction %in% c(
                "dnT",
                "gdT",
                "ILC",
                "MAIT"
            ) ~ "Unconventional T",
            wg2_scpred_prediction %in% c(
                "pDC",
                "cDC1",
                "cDC2",
                "ASDC"
            ) ~ "Dendritic",
            wg2_scpred_prediction %in% c(
                "CD14_Mono",
                "CD16_Mono"
            ) ~ "Monocyte",
            wg2_scpred_prediction %in% c(
                "HSPC"
            ) ~ "Other"
        ) %>% fct_relevel(c("B", "NK", "CD8 T", "CD4 T", "Unconventional T", "Dendritic", "Monocyte", "Other"))
    )


# ⚙️ ggplot themes ----

# define a ggplot theme that we can re-use across all the plots for consistent styling
# work in progress

theme_tenk10k <- function() {
    font <- "Helvetica"

    theme_classic() %+replace%
        theme(
            plot.title = element_blank(),
        )
}
