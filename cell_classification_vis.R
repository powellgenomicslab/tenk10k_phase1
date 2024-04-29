# get the pretrained model from scpred

# singularity shell /home/blabow/cell_classification.sif --bind /directflow
scpred_model <- readRDS("/hier_scpred.RDS")
saveRDS(scpred_model, "./scpred_tree.RDS")

# exit singularity

library(Seurat)
# library(DiagrammeR)
library(data.tree)
library(igraph)
library(tidyverse)
library(tidygraph)
library(ggraph)

scpred_model <- readRDS("/home/blabow/scpred_tree.RDS")

scpred_ggraph <- scpred_model %>%
  as_tbl_graph() %>%
  ggraph(layout = "igraph", algorithm = "tree") +
  geom_edge_diagonal() +
  geom_node_label(aes(label = name), size = 2) +
  coord_flip()
scpred_ggraph %>%
  ggsave(filename = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/scpred_tree.png", width = 10, height = 10)
