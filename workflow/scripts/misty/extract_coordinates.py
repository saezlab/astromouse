# %%
import scanpy as sc
import pandas as pd

# %%
if 'snakemake' in locals():
    tissue = snakemake.wilcards[0]
    adata_fp = snakemake.input[0]
    output_fp = snakemake.output[0]
else:
    tissue = 'brain'
    adata_fp = 'data/working/ST/{0}_wImages.h5ad'.format(tissue)
    output_fp = 'test.csv'


# %%
adata = sc.read_h5ad(adata_fp)

# %%
coords = adata.obs.filter(['array_row', 'array_col', 'library_id'])
coords[['x','y']] = pd.DataFrame(adata.obsm['spatial'], columns=['x', 'y'], index = coords.index)

coords.to_csv(output_fp)



# %%
