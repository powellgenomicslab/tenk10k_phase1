import sys
import numpy as np
import pandas as pd
import cna
import multianndata as mad

# celltype = sys.argv[1]
# resolution = sys.argv[2]
# analysis_name = sys.argv[3]

resolution = "major_cell_types"
analysis_name = "no_expr_pc_covars"
celltype = "Monocyte"

# read in the single cell data
madata = cna.read(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/multianndata/{resolution}/{analysis_name}/{celltype}_scDataObject.dimreduc.pca.gena_pheno.h5ad"
)

variant = "12:69350234:C:A"

# read in and parse the top enriched gene sets
fgsea = pd.read_csv(
    f"/directflow/SCCGGroupShare/projects/blabow/tenk10k_phase1/data_processing/csa_qtl/output/fgsea/{celltype}_npheno_{variant}_fgsea_sig.tsv",
    sep="\t",
)

# Parse the fgsea results turn gene set back to 
fgsea["leadingEdge"] = fgsea["leadingEdge"].str.split(",")
