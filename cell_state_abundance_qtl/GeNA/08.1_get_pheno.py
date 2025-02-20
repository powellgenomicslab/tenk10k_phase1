# This script uses methods described in https://github.com/immunogenomics/GeNA-applied/blob/main/run_gwas/define_phenotype.py
# Defines csaQTL sample-level and neighbourhood-level phenotypes, using CNA

import sys
import numpy as np
import pandas as pd
import cna
import multianndata as mad

celltype = sys.argv[1]
resolution = sys.argv[2]
analysis_name = sys.argv[3]
covs = sys.argv[4]

# resolution = "major_cell_types"
# analysis_name = "no_expr_pc_covars"
# celltype = "Dendritic"
# covs = "sex,age,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6,geno_PC7,BioHEART"

# input files
GeNA_dir = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/"
gena_sumstats_lead_snps_file = f"{GeNA_dir}/GeNA_sumstats_lead_snps_MAF_0.05.tsv"
genotypes_sig_file = f"{GeNA_dir}/GeNA_sig_snps.raw"
madata_file = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.pca.h5ad"

# parse arguments
covs_list = covs_list = covs.split(",")

# read in the lead SNPs
gena_sumstats_lead_snps = pd.read_csv(
    gena_sumstats_lead_snps_file, sep="\t"
).sort_values("P")

# read in genotypes for significant snps including non-lead snps
genotypes_sig = (
    pd.read_csv(genotypes_sig_file, sep="\t")
    .set_index("IID")
    .drop(["FID", "PAT", "MAT", "SEX", "PHENOTYPE"], axis=1)
    # .transpose()
)

# genotypes_sig.index = [id.split("_")[0] for id in genotypes_sig.index]
genotypes_sig.columns = [id.split("_")[0] for id in genotypes_sig.columns]

# filter genotypes to get just the lead snps, now col order is in order of significance
genotypes_lead_snps = genotypes_sig.loc[:, gena_sumstats_lead_snps["ID"]]

# read in the mult-anndata object for cell type
madata = cna.read(madata_file)

# add genotypes to sample metadata slot
madata.samplem = madata.samplem.join(genotypes_lead_snps, how="left")

# Store neighborhood-level phenotype for each SNP
npheno = pd.DataFrame({})

# copy uns
uns = madata.uns.copy()

# iterate through each lead SNP, calculate the sample-level and neighbourhood-level phenotype, corrected for covariates
for variant in gena_sumstats_lead_snps["ID"].tolist():
    # read in the value for k's (number of NAM-PC's used by GeNA)
    # NOTE: in actual GeNA model, k can be different for each SNP
    # Change this to read in the sumstats and use the actual correct k value in the future
    # variant = gena_sumstats_lead_snps["ID"][0]

    k = gena_sumstats_lead_snps.loc[
        gena_sumstats_lead_snps["ID"] == variant, "k"
    ].values[0]
    P = gena_sumstats_lead_snps.loc[
        gena_sumstats_lead_snps["ID"] == variant, "P"
    ].values[0]

    res = cna.tl.association(
        madata, madata.samplem[variant], covs=madata.samplem[covs_list], ks=[k]
    )

    # calculate the "sample-level phenotypes"
    spheno = pd.Series(
        np.dot(madata.uns["NAM_sampleXpc"].iloc[:, : res.k], res.beta),
        index=madata.uns["NAM_sampleXpc"].index,
        name=f"spheno_{variant}",
    )
    # store sample-level phenotypes in samplem
    madata.samplem = madata.samplem.join(spheno, how="left")
    npheno[f"npheno_{variant}"] = res.ncorrs

# add neighbourhood level phenotypes to the scanpy object
npheno.index = madata.obs.index
madata.obs = madata.obs.join(npheno, how="left")

# export sample and neighborhood level phenotypes with metadata
madata.samplem.to_csv(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_sample_pheno.csv"
)
madata.obs.to_csv(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_neighbourhood_pheno.csv"
)

# export the maultianndata object with phenotypes
madata.write(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_scDataObject.dimreduc.pca.gena_pheno.h5ad"
)
