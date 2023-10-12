# %%
import scanpy as sc
import pandas as pd

import matplotlib.pyplot as plt

# %%
import seaborn as sns

# %%
# Define input and output files
if 'snakemake' in locals():
    tissue = snakemake.wildcards[0]

    cellprop_cutoff = snakemake.params[0]

    adata_fp = snakemake.input[0]
    ctProp_fp = snakemake.input[1]
    ct_annot_fp = snakemake.input[2]

    output_fp = snakemake.output[0]

else:
    tissue = 'brain'
    
    sig_thresh = 0.05
    cellprop_cutoff = 0.1

    adata_fp = 'results/ST/{0}_wImages.h5ad'.format(tissue)
    ctProp_fp = 'results/ST/ST_{0}_deconvoluted.csv'.format(tissue)
    ct_annot_fp = 'data/original/MO/MO_cluster_metadata.csv'


# %%
#Load anndata
adata = sc.read_h5ad(adata_fp)
print(adata)

# %%
#Load celltype cluster annotation
annot = pd.read_csv(ct_annot_fp)
annot

# %%
#Load celltype proportionss from stereoscope
ctProp = pd.read_csv(ctProp_fp, index_col=0)
ctProp = ctProp.where(ctProp >= cellprop_cutoff) #set values below cutoff to NA
ctProp = pd.merge(ctProp, adata.obs.filter(['library_id'], axis = 1), left_index=True, right_index=True)

# %%
#Count number of of non-NA spots per celltype and sample
counts = ctProp.groupby(['library_id']).count().reset_index()
counts.columns = ['library_id'] + annot['clusterAbrv'].tolist()
counts = pd.melt(counts, id_vars='library_id', var_name='celltype', value_name='count')
counts = pd.merge(adata.obs.filter(['library_id', 'mouse', 'condition'], axis=1).drop_duplicates().reset_index(drop = True), counts, on='library_id')

# %%
#Median proportion of celltypes across non-NA spots per sample
medians = ctProp.groupby(['library_id']).median().reset_index()
medians.columns = ['library_id'] + annot['clusterAbrv'].tolist()
medians = pd.melt(medians, id_vars='library_id', var_name='celltype', value_name='median')
medians = pd.merge(adata.obs.filter(['library_id', 'mouse', 'condition'], axis=1).drop_duplicates().reset_index(drop = True), medians, on='library_id')

# %%
colors = sns.color_palette("tab10").as_hex()[0:6]

fig, axes = plt.subplots(1, 2 , figsize=(16, 15), dpi = 300)
axes = axes.flatten()

#Plot counts and medians per sample for each celltype
for df, name, ax in zip([counts, medians], ['count', 'median'], axes):

    sns.boxplot(data = df, x=name, y= 'celltype', hue = 'condition', color = 'lightgrey', order=annot['clusterAbrv'].sort_values(), fliersize=0, ax=ax)

    for mouse, color in zip(sorted(df['mouse'].unique()), colors):
        mouse_toplot = df[df['mouse'] == mouse]
        sns.stripplot(x=name, y='celltype', hue='condition', dodge=True, palette=[color] * 2, marker='o', data=mouse_toplot, order=annot['clusterAbrv'].sort_values(), ax=ax)

    handles, labels = ax.get_legend_handles_labels()
    handles = handles[0:2] + handles[2:-1:2]
    labels = labels[0:2] + sorted(counts['mouse'].unique())
    ax.legend(handles, labels)

axes[0].set_title('Number of spots with celltype')
axes[1].set_title('Median proportions in spots')
axes[0].set_ylabel('')
axes[1].set_ylabel('')

plt.suptitle(tissue.capitalize() + ' stereoscope results')
plt.tight_layout()

if 'snakemake' in locals(): plt.savefig(snakemake.output[0])
