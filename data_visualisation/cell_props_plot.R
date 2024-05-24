library(tidyverse)
library(ggsci)

cell_types <- read_csv("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/celltypes_by_cohort.csv")
head(cell_types)
dim(cell_types)

PropPlot_celltypes_cohort <- cell_types %>%
    ggplot(aes(x = cohort, fill = wg2_scpred_prediction)) +
    geom_bar(position = "fill") +
    theme_bw() +
    labs(y = "Cell type proportion (%)") +
    scale_fill_igv()

PropPlot_celltypes_cohort %>%
    ggsave(
        filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/scanpy/output/integrated_objects/figures/cell_types_proportions.png",
        width = 9, height = 5
    )
