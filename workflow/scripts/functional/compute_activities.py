# %%
import scanpy as sc
import pandas as pd
import decoupler as dc

# %%
if 'snakemake' in locals():
    tissue = snakemake.wildcards[0]
    network = snakemake.wildcards[1]
    adata_fp = snakemake.input[0]
    output_fp = snakemake.output[0]
    normalisation = snakemake.params[0]
    top_genes = int(snakemake.params[1])
    TF_conf = snakemake.params[2]
    method = snakemake.params[3]
else:
    tissue = 'brain'
    network = 'TFs'
    adata_fp = 'data/working/ST/{0}_wImages.h5ad'.format(tissue)
    output_fp = 'test.csv'
    normalisation = 'log1p'
    top_genes = int('300')
    TF_conf = 'ABC'
    method = 'mlm'


# %%
adata = sc.read_h5ad(adata_fp)

# %%
if normalisation == 'log1p':
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
elif normalisation == 'SCT':
    adata.layers['counts'] = adata.X
    adata.X = adata.layers['SCT'].copy()
    del adata.layers['SCT']
    adata

# %%
if network == 'pathways':
    model = dc.get_progeny(organism='mouse', top=top_genes)
elif network == 'TFs':
    model = dc.get_dorothea(organism='mouse', levels=[c for c in TF_conf])
else:
    raise ValueError('The "network" wildcard can only take on "pathways" or "TFs" as value, to run either Progeny or Dorothea regulons')

# %%
dc.decouple(mat=adata, net=model, source='source', target='target', weight='weight', methods = method,  verbose=True, use_raw=False)

# %%
print(adata)

# %%
acts = dc.get_acts(adata, obsm_key='{0}_estimate'.format(method))
print(acts)

# %%
acts = pd.DataFrame(acts.X, index=acts.obs.index, columns=acts.var.index)

# %%
acts.head()

# %%
acts.to_csv(output_fp)


