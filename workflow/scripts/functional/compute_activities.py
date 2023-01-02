# %%
import scanpy as sc
import pandas as pd
import decoupler as dc
import os

# %%
if 'snakemake' in locals():
    tissue = snakemake.wildcards[0]
    net_fp = snakemake.input.get('net', '')

    network = snakemake.wildcards[1]
    adata_fp = snakemake.input[0]
    output_fp = snakemake.output[0]
    
    conf = snakemake.params[0]

else:
    tissue = 'brain'
    network = 'TFs'
    adata_fp = 'data/working/ST/{0}_wImages.h5ad'.format(tissue)
    output_fp = 'test.csv'
    
    conf = {'normalisation': 'log1p', 'top_targets': 300, 'method': 'mlm'}


# %%
adata = sc.read_h5ad(adata_fp)

# %%
#Defining normalisation
if conf.get('normalisation') == 'log1p':
    print('INFO: using log1p normalised counts')
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
elif conf.get('normalisation') == 'SCT':
    print('INFO: using SCT normalised counts')
    adata.layers['counts'] = adata.X
    adata.X = adata.layers['SCT'].copy()
    del adata.layers['SCT']
elif conf.get('normalisation') is None:
    raise ValueError('The config file is missing a normalisation method for {0}. Set it to either \'log1p\' or \'SCT\'.'.format(network))
else:
    raise ValueError('The normalisation method \'{0}\' is not implemented. Set it to either \'log1p\' or \'SCT\'.'.format(conf.get('normalisation')))

# %%
# Load regulon network
if network == 'pathways':
    model = dc.get_progeny(organism='mouse', top=conf.get('top_targets'))
elif network == 'TFs':
    model = dc.get_dorothea(organism='mouse', levels=[c for c in conf.get('levels')])
elif network == 'GRNs':
    if net_fp == '':
        raise ValueError('No file was provided for the GRN regulons!')
    else:
        model = pd.read_csv(os.path.join(net_fp, tissue + '.csv'), sep=',')
else:
    raise ValueError('The "network" wildcard can only take on "pathways", "GRNS" or "TFs" as value, to run either Progeny, celloracle GRNs or Dorothea regulons')

# %%
#run decoupler on the given PKN
dc.decouple(mat=adata, net=model, source='source', target='target', weight='weight', methods = conf.get('method'),  verbose=True, use_raw=False)

# %%
print(adata)

# %%
#extract activities from the anndata object
acts = dc.get_acts(adata, obsm_key='{0}_estimate'.format(conf.get('method')))
print(acts)
acts = pd.DataFrame(acts.X, index=acts.obs.index, columns=acts.var.index)

# %%
acts.head()

# %%
acts.to_csv(output_fp)


