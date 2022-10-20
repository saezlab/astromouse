import pandas as pd
import numpy as np
import celloracle as co
import mudata
import scanpy as sc
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--path_mdata', required=True)
parser.add_argument('-b','--base_grn', required=True)
parser.add_argument('-l','--path_links', required=True)
args = vars(parser.parse_args())

path_mdata = args['path_mdata']
path_base_grn = args['base_grn']
path_links = args['path_links']
n_comps = 50

# Read rna adata
adata = mudata.read(os.path.join(path_mdata, 'rna'))

# Load base grn
base_GRN = pd.read_csv(path_base_grn, index_col=0)

# Instantiate Oracle object
oracle = co.Oracle()

# Find HVG
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, batch_key='bio_origin')
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
num_hvg_genes = 3000
batch_msk = np.array(adata.var.highly_variable_nbatches > 1)
hvg = adata.var[batch_msk].sort_values('highly_variable_nbatches').tail(num_hvg_genes).index
adata.var['highly_variable'] = [g in hvg for g in adata.var.index]
adata.var = adata.var[['highly_variable','highly_variable_nbatches']]
adata = adata[:,hvg]
adata.X = adata.layers['counts'].copy()
adata.obs['tissue'] = 'tissue'

# Instantiate Oracle object
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="tissue",
                                   embedding_name="X_umap")

# You can load TF info dataframe with the following code
oracle.import_TF_data(TF_info_matrix=base_GRN)

# Perform PCA
oracle.perform_PCA()

# KNN imputation
n_cell = oracle.adata.shape[0]
print("cell number is: {0}".format(n_cell))

k = int(0.025 * n_cell)
print('Auto-selected k is: {0}'.format(k))

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=os.cpu_count())

# Compute GRN
links = oracle.get_links(cluster_name_for_GRN_unit="tissue", alpha=10,
                         verbose_level=10)

# Save Links object
links.to_hdf5(file_path=path_links)

