# mm activate gena-env

import sys
import seaborn as sns
import pandas as pd
import numpy as np
import multianndata as mad
import cna
import scanpy as sc
from matplotlib import pyplot as plt

sys.path.append(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/cell_state_abundance_qtl/GeNA"
)
import pgen  # pgen reader functions from tensorQTL

np.random.seed(0)

# celltype = sys.argv[1]
# resolution = sys.argv[2]

resolution = "major_cell_types"
analysis_name = "no_expr_pc_covars"
covs = "sex,age,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6,geno_PC7,BioHEART"

# Dendritic cells
# celltype = "Dendritic"
# variant_ids = "12:121634368:T:G,22:19725903:TCG:T,3:195961658:A:ACAGGGGGCGGAGCTCACAGTAATGCTCTGG"  # lead variant intersecting immune deficiency genes

# # Monocytes
# celltype = "Monocyte"
# variant_ids = "12:121634368:T:G,22:19725903:TCG:T,16:30169425:T:TG"  # lead variant intersecting immune deficiency genes

# Unconventional_T
# celltype = "Unconventional_T"
# variant_ids = "12:121634284:G:T"  # lead variant intersecting immune deficiency genes

# NK cells
# celltype = "NK"
# variant_ids = "19:16331208:G:A"  # lead variant intersecting T cell Abundance GWAS

# B cells
celltype = "B"
variant_ids = "14:105840505:T:C"

## REAL
plink_prefix_path = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15"
## PERMUTED
# plink_prefix_path = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted/merged_common_variants_standard_chr_geno_0.15_singlecell_cohort_perm1"

# parse args
variant_id_list = variant_ids.split(",")
covs_list = covs.split(",")

# output directories
outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)
figdir = f"{outdir}/figures"


# variants passing genome-wide significance threshold from the REAL test
gena_sumstats_sig_path = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/GeNA_sumstats_sig_MAF_0.05.csv"
gena_sumstats_sig = pd.read_csv(gena_sumstats_sig_path)
gena_sumstats_sig["permuted"] = "Real genotype csaQTL"
# variants passing genome-wide significance threshold from the PERMUTED test
# gena_sumstats_sig_perm_path = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}_perm/GeNA_sumstats_sig.csv"
# gena_sumstats_sig_perm = pd.read_csv(gena_sumstats_sig_perm_path)
# gena_sumstats_sig_perm["permuted"] = "Permuted genotype csaQTL"

# gena_sumstats_sig = gena_sumstats_sig_perm

# check how many intersecting between permuted and real genotype data
# len(set(gena_sumstats_sig["ID"]).intersection(set(gena_sumstats_sig_perm["ID"])))

madata = cna.read(
    f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.pca.h5ad"
)

# top 10 significant variants
variant_df = gena_sumstats_sig.sort_values(by=["P"]).loc[
    gena_sumstats_sig["ID"].isin(variant_id_list), ["ID", "k", "P", "permuted"]
]

# read in the genotypes and add them to madata object
pgr = pgen.PgenReader(plink_prefix_path)
dosages = pgr.read_dosages_list(variant_df["ID"].tolist()).transpose()


# ################################
# dosages.index.duplicated().sum()
# # check that samples in plink file match samples in madata
# len(set(madata.samplem.index).intersection(set(dosages.index)))
# len(set(madata.samplem.index))
# testdf = madata.samplem.join(dosages, how="left")

# len(
#     set(madata.samplem["sample_perm0"]).difference(set(madata.samplem["sample_perm1"]))
# )
# len(set(madata.samplem["sample_perm2"]).difference(set(madata.samplem["sample_perm3"])))

# ######################

madata.samplem = madata.samplem.join(dosages, how="left")

for variant in variant_df["ID"].tolist():
    # read in the value for k's (number of NAM-PC's used by GeNA)
    # NOTE: in actual GeNA model, k can be different for each SNP
    # Change this to read in the sumstats and use the actual correct k value in the future
    # ks = pd.read_csv(
    #     f"{outdir}/output/GeNA/{resolution}/{celltype}/{analysis_name}/ks.csv",
    #     names=["k"],
    # )
    # kmax = ks.max()[0]

    k = variant_df.loc[variant_df["ID"] == variant, "k"].values[0]
    P = variant_df.loc[variant_df["ID"] == variant, "P"].values[0]
    permuted = variant_df.loc[variant_df["ID"] == variant, "permuted"].values[0]

    res = cna.tl.association(
        madata, madata.samplem[variant], covs=madata.samplem[covs_list], ks=[k]
    )
    # calculate the "sample-level phenotypes"
    madata.samplem["Sample-level phenotype"] = np.dot(
        madata.uns["NAM_sampleXpc"].iloc[:, : res.k], res.beta
    )
    sample_pheno = madata.samplem[["Sample-level phenotype", variant]]
    sample_pheno.index = madata.uns["NAM_sampleXpc"].index

    # Convert dosage from numeric encoding to actual genotypes
    # TODO Maybe better to use genotype rather than dosage? not sure if dosage is always E {0,1,2}
    # genotypes_dict = {v: v.split(":")[2:4] for v in variant_id_list}
    # some dosage are -9, which is for missing

    # def geno_to_letters(ref, alt, geno):
    #     if geno == 0:
    #         out = f"{ref}/{ref}"
    #     elif geno == 1:
    #         out = f"{ref}/{alt}"
    #     elif geno == 2:
    #         out = f"{alt}/{alt}"
    #     elif geno == -9:
    #         out = "Missing genotype"
    #     return out

    # ref, alt = genotypes_dict[variant]

    # sample_pheno["Genotype"] = [
    #     geno_to_letters(ref, alt, g) for g in sample_pheno[variant]
    # ]

    # plot_order = [f"{ref}/{ref}", f"{ref}/{alt}", f"{alt}/{alt}"]

    sns.set_theme(style="darkgrid")
    ax = sns.boxplot(
        data=sample_pheno,
        x=f"{variant}",
        y="Sample-level phenotype",
        showfliers=False,
        # order=plot_order
    )
    ax = sns.swarmplot(
        data=sample_pheno,
        x=f"{variant}",
        y="Sample-level phenotype",
        color="grey",
        size=1,
        # order=plot_order,
    ).set(title=f"{celltype}, P = {P}, {permuted}")
    plt.savefig(
        fname=f"{figdir}/{resolution}/boxplots/intersections_real_geno/{celltype}_{analysis_name}_{variant}_intersecting_curated_deficiency_genes_box_plot.png",
        dpi=300,
    )
    plt.clf()

# Todo: order plot by homoref, het, homoalt
#


# ---------------------------------------------------------------------------
# convert genotypes for plotting
# genotypes = pgr.read_list(variant_id_list).transpose()

# # convert genotypes from numeric encoding
# genotypes_letters = pd.DataFrame(index=genotypes.index)
# for v in variant_id_list:
#     ref, alt = genotypes_dict[v]
#     genotypes_letters[v] = [geno_to_letters(ref, alt, g) for g in genotypes[v]]
