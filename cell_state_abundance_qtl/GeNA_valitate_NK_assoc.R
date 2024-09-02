library(data.table)

sumstats <- fread("/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/GeNA_sumstats.txt")

# intersect the GeNA significant associations ? or lead SNPs
sumstats[, ID := paste0(sumstats$`#CHROM`, ":", sumstats$POS)] ## Fix the ID column - I think the ID must be old pos from GrCH37 ?? 

