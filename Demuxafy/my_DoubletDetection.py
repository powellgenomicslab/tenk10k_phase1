import os
import numpy as np
import doubletdetection
import tarfile
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import sys
import pandas as pd

# Load read10x function from mods directory

mods_path = "/opt/Demultiplexing_Doublet_Detecting_Docs/mods" ## Do not change - this is the path to the mods folder in the singularity image with custom script for loading 10x data in python
sys.path.append(mods_path)
import read10x

### Set up parameters and variables ###
cellranger_dir = str(sys.argv[1])
outdir = str(sys.argv[2])

if not os.path.isdir(outdir):
  os.mkdir(outdir)


### Read in data ###
raw_counts = read10x.import_cellranger_mtx(cellranger_dir + "/matrix.mtx.gz")

try:
  barcodes_df = read10x.read_barcodes(cellranger_dir + "/barcodes.tsv.gz")
except:
  try:
    barcodes_df = read10x.read_barcodes(cellranger_dir + "/barcodes.tsv")
  except:
    print("No barcode file in provided counts matrix directory. Please double check the directory or provide the full path to the barcode file to use.")

print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Remove columns with all 0s
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]
print('Counts matrix shape after removing unexpressed genes: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

clf = doubletdetection.BoostClassifier(n_iters=50, clustering_algorithm='phenograph', standard_scaling=True, verbose = True)
doublets = clf.fit(raw_counts).predict(p_thresh=1e-16, voter_thresh=50)

results = pd.Series(doublets, name="DoubletDetection_DropletType")
dataframe = pd.concat([barcodes_df, results], axis=1)
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")

dataframe.to_csv(os.path.join(outdir,'DoubletDetection_doublets_singlets.tsv'), sep = "\t", index = False)


### Figures ###
doubletdetection.plot.convergence(clf, save=os.path.join(outdir,'convergence_test.pdf'), show=False, p_thresh=1e-16, voter_thresh=0.5)

f3 = doubletdetection.plot.threshold(clf, save=os.path.join(outdir,'threshold_test.pdf'), show=False, p_step=6)


### Make summary of singlets and doublets and write to file ###
summary = pd.DataFrame(dataframe.DoubletDetection_DropletType.value_counts())
summary.index.name = 'Classification'
summary.reset_index(inplace=True)
summary = summary.rename({'DoubletDetection_DropletType': 'Droplet N'}, axis=1)

summary.to_csv(os.path.join(outdir,'DoubletDetection_summary.tsv'), sep = "\t", index = False)
