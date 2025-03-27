# mm activate gena-env

import sys
import pandas as pd
import numpy as np

np.random.seed(0)

celltype = sys.argv[1]
resolution = sys.argv[2]
analysis_name = sys.argv[3]

# resolution = "major_cell_types"
# analysis_name = "no_expr_pc_covars"
# celltype = "Dendritic"

GeNA_dir = f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/GeNA/{resolution}/{celltype}/{analysis_name}/"

# plink_prefix_path = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15"
plink_prefix_path = f"{GeNA_dir}/GeNA_sig_snps"
# output directory
output = f"{GeNA_dir}/GeNA_sumstats_lead_snps_MAF_0.05.tsv"


# variants passing genome-wide significance threshold from the REAL test
gena_sumstats_sig_path = f"{GeNA_dir}/GeNA_sumstats_sig_MAF_0.05.csv"
gena_sumstats_sig = pd.read_csv(gena_sumstats_sig_path)  # .set_index("ID")
gena_sumstats_sig = gena_sumstats_sig.sort_values("P")

# read in genotypes, transform it to have 1 snp per row and one individual per column
genotypes_sig = (
    pd.read_csv(f"{plink_prefix_path}.raw", sep="\t")
    .set_index("IID")
    .drop(["FID", "PAT", "MAT", "SEX", "PHENOTYPE"], axis=1)
    .transpose()
)
genotypes_sig.index = [
    id.split("_")[0] for id in genotypes_sig.index
]  # remove _refallele encoding from ids

# sort genotypes df in order of significance
genotypes_sig = genotypes_sig.reindex(gena_sumstats_sig["ID"])

# use pruning logic from GeNA paper
# Reference: https://github.com/immunogenomics/GeNA-applied/blob/main/run_gwas/define_loci.py

res = gena_sumstats_sig
G = genotypes_sig
all_loci = pd.DataFrame({})
window_sz = 1000000  # 1 MB

while res.shape[0] > 0:
    new = res.iloc[0:1, :]
    new.insert(0, "celltype", celltype)
    all_loci = pd.concat([all_loci, new])  # Add lead snp for new locus
    lead_POS = res.POS.values[0]
    res = res.drop(index=res.index[0])  # Remove lead snp
    if res.shape[0] > 0:  # check if any snps left
        G_corrs = (
            np.abs(G.T.corr()).iloc[1:, :] ** 2
        )  # Compute R^2 to other remaining snps
        G = G.drop(index=G.index[0])  # Remove lead snp
        i_drop = np.where(G_corrs.iloc[:, 0].values > 0.8)[
            0
        ]  # Find any other snps in locus
        res = res.drop(index=res.index[i_drop])  # Rm other snps in locus
        G = G.drop(index=G.index[i_drop])  # Rm other snps in locus

        # Also rm anything in 1MB window around lead snp
        in_window = (res.POS > lead_POS - (window_sz / 2)) & (
            res.POS < lead_POS + (window_sz / 2)
        )
        i_drop = np.where(in_window)[0]
        res = res.drop(index=res.index[i_drop])  # Rm other snps in locus
        G = G.drop(index=G.index[i_drop])  # Rm other snps in locus

all_loci.reset_index()
all_loci.to_csv(output, sep="\t", index=False)
