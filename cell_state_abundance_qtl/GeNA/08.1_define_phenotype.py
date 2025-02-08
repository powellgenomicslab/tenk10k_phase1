# Script taken from from https://github.com/immunogenomics/GeNA-applied/blob/adaeb11e1eef070e0b2d9b70ba94458cb8a5c5a6/run_gwas/define_phenotype.py
# This script defines the csaQTL phenotype associated with each SNP in a provided
# set. Specifically, CNA is applied to compute and store per SNP the associated
# phenotype values per sample and across neighborhoods. Additionally, the correlation
# between expression per variable gene and neighborhood phenotype is stored.
#
# Inputs
#     sc_object: filepath to a single-cell data object in multianndata format
#           Expected to include alternative allele genotype doses for all snps of interest in the samplem dataframe
#     covs: covariates included in the GeNA GWAS, formatted as "cov1,cov2,cov3"; must correspond to columns in the samplem dataframe within the single-cell data object
#     loci: filepath to a GeNA summary statistics file; phenotypes will be defined for all SNPs in the file (usually these are previously refined to represent lead SNPs)
#     outfile_path: filepath to a location where all results objects will be stored
#
# Outputs
#     spheno.tsv: a text file with one row per individual in the cohort and one column per SNP storing sample-level phenotype values
#     npheno.tsv: a text file with one row per neighborhood and one column per SNP storing neighborhood-level phenotype values
#         Neighborhoods are labeled by the ID of their anchor cell. If batch was included as a covariate, batchy neighborhoods are dropped.
#     cna_res_[SNP].p: a pickled CNA results object (details at https://github.com/immunogenomics/CNA) per SNP
#         Additionally contains 'vargene_cors' attribute, a dataframe storing the Pearson's correlation to neighborhood phenotype per variable gene

import argparse
import pandas as pd
import numpy as np
import cna, pickle
import statsmodels.api as sm

np.random.seed(0)

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--loci", type=str)
parser.add_argument("--sc_object", type=str)
parser.add_argument("--outfile_path", type=str)
parser.add_argument("--covs", type=str, default=None)
args = parser.parse_args()
print("\n\n****")
print(args)
print("****\n\n")

all_loci = pd.read_table(args.loci)
d = cna.read(args.sc_object)

if args.covs is not None:
    args.covs = args.covs.split(",")
    args.covs = d.samplem[args.covs]

# Store sample-level phenotype for each SNP
spheno = pd.DataFrame({})
spheno.insert(0, "id", d.samplem.index)

# Store neighborhood-level phenotype for each SNP
npheno = pd.DataFrame({})

# Iterate through SNPs
for i_allele in np.arange(all_loci.shape[0]):
    sel_k = all_loci.k.values[i_allele]  # model selected by GeNA
    G_allele = all_loci.ID.values[i_allele]
    res = cna.tl.association(
        d, d.samplem[G_allele], covs=args.covs, batches=d.samplem.batch, ks=[sel_k]
    )

    # Variable gene correlations
    vargene_cors = []
    for i_gene in np.arange(d.var.shape[0]):
        vargene_cors.append(np.corrcoef(d.X[res.kept, i_gene], res.ncorrs)[0, 1])
    res.vargene_cors = pd.DataFrame({"gene": d.var.index, "cor": vargene_cors})

    # Save CNA results object
    pickle.dump(res, open(args.outfile_path + "cna_res_" + G_allele + ".p", "wb"))

    # Save phenotype values per neighborhood
    npheno[G_allele] = res.ncorrs

    # Define and save phenotype values per sample
    d.samplem["pheno"] = np.dot(d.uns["NAM_sampleXpc"].iloc[:, : res.k], res.beta)
    spheno[G_allele] = d.samplem.pheno.values

npheno.insert(0, "cellid", d.obs.index[res.kept])
npheno.to_csv(args.outfile_path + "npheno.tsv", index=False, sep="\t")
spheno.to_csv(args.outfile_path + "spheno.tsv", index=False, sep="\t")
