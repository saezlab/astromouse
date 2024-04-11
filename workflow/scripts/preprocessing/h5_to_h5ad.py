# %%
from scipy.sparse import csr_matrix
import h5py
import pandas as pd

import platform
import os

import scanpy as sc
import numpy as np

# %%
if 'snakemake' in locals():
    input_file = snakemake.input[0]
    wildcards = snakemake.wildcards
    
    main_assay = 'SCT'

else:
    input_file = 'test.h5'
    main_assay = 'SCT'

# %%
def readh5_to_anndata(file, main_assay = 'RNA'):
    """ Reads from h5 file and returns anndata object

    Currently only support a single type of data (i.e. data OR counts OR scaled data from Seurat).
    Args:
        input_file (path): Path to h5 file
        main_assay (string): Name of assay in h5 file to use as X

    Returns:
        Anndata: Anndata with main_assay as X and all other assays as layers.
    """
    with h5py.File(file, "r") as f:

        # check that main assay is in f
        if main_assay not in list(f['assays'].keys()):
            raise ValueError('Main assay %s not in h5 file' % main_assay)
        
        # check that obs is in f
        if 'obs' not in list(f.keys()):
            raise ValueError('Obs not in h5 file')

        # create anndata object
        print('Loading main assay %s' % main_assay)
        X = f['assays'][main_assay]

        mat = csr_matrix((X['values'], X['indices'], X['indptr']), shape=(X['dims'][0], X['dims'][1]))
        adata = sc.AnnData(X=mat, dtype=np.float64)

        print('INFO: first 20 cells and genes in the X matrix')
        print(adata.X[0:20, 0:20])

        adata.obs_names = X['obs_names'][()].astype(str)
        adata.var_names = X['var_names'][()].astype(str)

        # add obs
        print('Loading obs')
        obs = f['obs']
        for key in obs.keys():
            # Check if key is a factor and if so convert
            if key + '_levels' in obs.keys():
                idx = obs[key][()]
                adata.obs[key] = obs[key + '_levels'][()][idx].astype(str)#.str.decode("utf-8")
            elif '_levels' not in key and 'colnames' not in key:
                # if bytes in [x for x in obs[key][()].apply(type).unique()]:]
                #     adata.obs[key] = obs[key][()].astype(str)#.str.decode("utf-8")
                # else:
                adata.obs[key] = obs[key][()]

        print('INFO: head of obs')
        print(adata.obs.head())

        # load remaining assays as layers
        layers = [key for key in f['assays'].keys() if main_assay not in key]
        adata.layer = {}
        for layer in layers:
            print('Loading layer: %s' % layer)
            X = f['assays'][layer]

            # add only if dims match
            if X['dims'][0] == adata.shape[0] and X['dims'][1] == adata.shape[1]:
                adata.layer[layer] = csr_matrix((X['values'], X['indices'], X['indptr']), shape=(X['dims'][0], X['dims'][1]))
            else:
                raise Warning('Layer %s not added because dims do not match' % layer)
            print('INFO: first 20 cells and genes in the %s layer' % layer)
            print(adata.layer[layer][0:20, 0:20])

        # load reductions
        if 'reductions' in f.keys():
            adata.obsm = {}
            print('Loading reductions')
            reductions = f['reductions']
            for red in reductions.keys():
                print('\tLoading reduction: %s' % red)
                embedding = reductions[red]

                adata.obsm['X_' + red] = np.array([embedding[key][()] for key in embedding.keys()], dtype=np.float64).T

    return adata

# %%
adata = readh5_to_anndata(input_file, main_assay)

if 'snakemake' in locals():
    adata.write_h5ad(snakemake.output[0])
