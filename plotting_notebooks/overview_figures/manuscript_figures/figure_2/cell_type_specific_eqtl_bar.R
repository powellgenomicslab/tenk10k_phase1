library(tidyverse)
library(glue)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(qvalue)

source("plotting_notebooks/overview_figures/manuscript_figures/tenk_data_vis_utils.R")

# significance threshold to use
qval_threshold <- 0.05

# ⚙️ functions ----

# Code adpated from the STAR package https://github.com/xihaoli/STAAR/blob/dc4f7e509f4fa2fb8594de48662bbd06a163108c/R/CCT.R wtih a modifitcaiton: when indiviudal p-value = 1, use minimum p-value
#' An analytical p-value combination method using the Cauchy distribution
#'
#' The \code{CCT} function takes in a numeric vector of p-values, a numeric
#' vector of non-negative weights, and return the aggregated p-value using Cauchy method.
#' @param pvals a numeric vector of p-values, where each of the element is
#' between 0 to 1, to be combined.
#' @param weights a numeric vector of non-negative weights. If \code{NULL}, the
#' equal weights are assumed.
#' @return the aggregated p-value combining p-values from the vector \code{pvals}.
#' @examples pvalues <- c(2e-02, 4e-04, 0.2, 0.1, 0.8)
#' @examples CCT(pvals = pvalues)
#' @references Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association 115}(529), 393-402.
#' (\href{https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1554485}{pub})
#' @export

CCT <- function(pvals, weights = NULL) {
    #### check if there is NA
    if (sum(is.na(pvals)) > 0) {
        stop("Cannot have NAs in the p-values!")
    }

    #### check if all p-values are between 0 and 1
    if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
        stop("All p-values must be between 0 and 1!")
    }

    #### check if there are p-values that are either exactly 0 or 1.
    is.zero <- (sum(pvals == 0) >= 1)
    is.one <- (sum(pvals == 1) >= 1)
    # if(is.zero && is.one){
    #  stop("Cannot have both 0 and 1 p-values!")
    # }
    if (is.zero) {
        return(0)
    }
    if (is.one) {
        # warning("There are p-values that are exactly 1!")
        return(min(1, (min(pvals)) * (length(pvals))))
    }

    #### check the validity of weights (default: equal weights) and standardize them.
    if (is.null(weights)) {
        weights <- rep(1 / length(pvals), length(pvals))
    } else if (length(weights) != length(pvals)) {
        stop("The length of weights should be the same as that of the p-values!")
    } else if (sum(weights < 0) > 0) {
        stop("All the weights must be positive!")
    } else {
        weights <- weights / sum(weights)
    }

    #### check if there are very small non-zero p-values
    is.small <- (pvals < 1e-16)
    if (sum(is.small) == 0) {
        cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
    } else {
        cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
        cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
    }

    #### check if the test statistic is very large.
    if (cct.stat > 1e+15) {
        pval <- (1 / cct.stat) / pi
    } else {
        pval <- 1 - pcauchy(cct.stat)
    }
    return(pval)
}

get_CCT_pvalue <- function(pvalue, weights = NULL) {
    pvals <- pvalue
    notna <- which(!is.na(pvals))
    if (length(notna) > 0) {
        pvals <- pvals[!is.na(pvals)]
        cctpval <- CCT(pvals, weights = weights)
    } else {
        cctpval <- NA
    }
    return(cctpval)
}


compute_fdr_celltype <- function(ct, eqtl_dir, save = FALSE) {
    out_file <- paste0(eqtl_dir, "acat_fdr_gene_files/", ct, "with_snp_pos.tsv")
    if (!file.exists(out_file)) {
        raw_pvals_file <- paste0(eqtl_dir, ct, "_gene_tests_raw_pvals.txt")
        library(data.table)

        # Read the TSV file into a data.table
        raw_pvals_dt <- fread(raw_pvals_file, header = FALSE)

        # Set the column names
        setnames(raw_pvals_dt, c("chrom", "snp_pos", "gene", "pval"))

        # Add new columns
        raw_pvals_dt <- raw_pvals_dt %>%
            mutate(pval_cct = get_CCT_pvalue(pval))
        raw_pvals_dt$qv <- qvalue(raw_pvals_dt$pval_cct, pi0 = 1)$qvalue
        return(raw_pvals_dt) # WRONG
    }
}


# compute FDR results ----

eqtl_dir <- "/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/associatr/"

celltypes <- get_latest_metadata() %>%
    pull(wg2_scpred_prediction) %>%
    unique()

qv_df <- compute_fdr_celltype(celltypes[1], eqtl_dir = eqtl_dir)


# 📚 read in gene-level ACAT + FDR corrected results ----




# compute FDR for the raw pvalues




#' read and combine fdr results for all the cells
# get_fdr_results <- function(celltype, qval_threshold) {

#     read_csv(glue("/directflow/SCCGGroupShare/projects/anncuo/TenK10K_pilot/tenk10k/eqtl_results/associatr/acat_fdr_gene_files/{celltype}.tsv")) %>%
#         filter(qv < qval_threshold) %>%
#         mutate(celltype = celltype) %>%
#         return()
# }

celltype_fdr <- celltypes %>%
    map(\(celltype) get_fdr_results(celltype = celltype, qval_threshold = qval_threshold)) %>%
    list_rbind()
