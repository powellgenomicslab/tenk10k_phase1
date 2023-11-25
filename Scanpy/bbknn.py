import bbknn
import numpy as np
import scanpy as sc

# load combined AnnData object (first 64 libraries)
output_dir = "/share/ScratchGeneral/anncuo/tenk10k/data_processing/celltypist/"
filename_all = output_dir+"all_samples_concatenated_celltypist_predicted.h5"
merge = sc.read(filename_all)

print(f'Total number of cells: '{merge.shape[0]})
print(f'Total number of genes: '{merge.shape[1]})
print(f'Number of samples to integrate: '{len(merge.obs['sample].unique())})

bbknn.bbknn(merge, batch_key='sample')
