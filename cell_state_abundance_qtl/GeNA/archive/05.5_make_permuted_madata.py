import sys
import pandas as pd
import multianndata as mad
import cna
import scanpy as sc
from matplotlib import pyplot as plt

# celltype = sys.argv[1]
# resolution = sys.argv[2]

celltype = "Dendritic"
resolution = "major_cell_types"

outdir = (
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl"
)

# read in the sample metadata (contains shuffled sample id mappings used in eQTL analysis)
sample_meta = pd.read_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/saige-qtl_tenk10k-genome-2-3-eur_input_files_241210_covariates_sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv"
)

# construct dict mapping between real ids and shuffled ids
perm0_map_dict = {
    sample_meta.loc[index, "sample_id"]: sample_meta.loc[index, "sample_perm0"]
    for index in sample_meta.index
}

madata = cna.read(
    f"{outdir}/data/h5/{resolution}/{celltype}_scDataObject.dimreduc.pca.h5ad"
)

# change sampleids in obs
madata.obs.id = [perm0_map_dict[id] for id in madata.obs["id"]]

# change sampleids in samplem
madata.samplem.index = [perm0_map_dict[id] for id in madata.samplem.index]


# psam_perm["#IID"] = [perm0_map_dict[id] for id in psam["#IID"]]
### graveyard ####


# attempt to rename plink .psam file with permuted sample ids
# should work but first I would need to filter the plink file to just include samples within the single cell cohort


import sys
import os
import pandas as pd

outdir = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/permuted"
in_plink = "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/plink/merged_common_variants_standard_chr_geno_0.15"

# -----------
#  Create permuted genotype files
# -----------

# read in the sample metadata (contains shuffled sample id mappings used in eQTL analysis)
sample_meta = pd.read_csv(
    "/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/data/saige-qtl_tenk10k-genome-2-3-eur_input_files_241210_covariates_sex_age_geno_pcs_shuffled_ids_tob_bioheart.csv"
)

# construct dict mapping between real ids and shuffled ids
perm0_map_dict = {
    sample_meta.loc[index, "sample_id"]: sample_meta.loc[index, "sample_perm0"]
    for index in sample_meta.index
}

if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

# read in the input psam file
psam = pd.read_csv(f"{in_plink}.psam", delimiter="\t")
psam_perm = psam.copy()

# replace real sample ids with the permuted ones
psam_perm["#IID"] = [perm0_map_dict[id] for id in psam["#IID"] if id in psam["#IID"]]

for id in psam["#IID"]:
    print(perm0_map_dict[id])
