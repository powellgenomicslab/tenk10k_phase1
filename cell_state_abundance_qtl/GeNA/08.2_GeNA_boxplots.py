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

# celltype = sys.argv[1]
# resolution = sys.argv[2]

celltype = "Dendritic"
resolution = "major_cell_types"
analysis_name = "no_expr_pc_covars"
## REAL
# plink_prefix_path = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15"
## PERMUTED
plink_prefix_path = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted/merged_common_variants_standard_chr_geno_0.15_singlecell_cohort_perm0"

# variant_ids = "5:24615560:T:G,2:121911769:T:TAA" # NK
# variant_ids = "1:143272081:G:C,13:73835066:C:CT"  # Dendritic
covs = "sex,age,geno_PC1,geno_PC2,geno_PC3,geno_PC4,geno_PC5,geno_PC6,geno_PC7,BioHEART"
variant_ids = "1:143272081:G:C,10:42301537:A:G,12:121634284:G:T,19:27243890:C:T"  # Dendritic no exprPCs real geno
# variant_ids = ""

# parse args
variant_id_list = variant_ids.split(",")
covs_list = covs.split(",")

outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)
figdir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/figures"

madata = cna.read(
    f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.pca.h5ad"
)


# read in the genotypes and add them to madata object
pgr = pgen.PgenReader(plink_prefix_path)
dosages = pgr.read_dosages_list(variant_id_list).transpose()


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

for variant in variant_id_list:
    # read in the value for k's (number of NAM-PC's used by GeNA)
    # NOTE: in actual GeNA model, k can be different for each SNP
    # Change this to read in the sumstats and use the actual correct k value in the future
    ks = pd.read_csv(
        f"{outdir}/output/GeNA/{resolution}/{celltype}/{analysis_name}/ks.csv",
        names=["k"],
    )
    kmax = ks.max()[0]

    res = cna.tl.association(
        madata, madata.samplem[variant], covs=madata.samplem[covs_list], ks=[kmax]
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

    # def geno_to_letters(ref, alt, geno):
    #     if geno == 0:
    #         out = f"{ref}/{ref}"
    #     elif geno == 1:
    #         out = f"{ref}/{alt}"
    #     elif geno == 2:
    #         out = f"{alt}/{alt}"
    #     return out

    # ref, alt = genotypes_dict[variant]

    # sample_pheno["Genotype"] = [
    #     geno_to_letters(ref, alt, g) for g in sample_pheno[variant]
    # ]

    # plot_order = [f"{ref}/{ref}", f"{ref}/{alt}", f"{alt}/{alt}"]
    # plot_order = [0, 1, 2]

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
        y="pheno",
        color="grey",
        size=1,
        # order=plot_order,
    )
    plt.savefig(
        fname=f"{figdir}/{resolution}/boxplots/{celltype}_{analysis_name}_{variant}_box_plot.png",
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
