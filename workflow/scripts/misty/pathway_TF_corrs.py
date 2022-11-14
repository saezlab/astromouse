# %%
import scanpy as sc
import pandas as pd
import decoupler as dc
import numpy as np
import os
import re

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# %%
from scipy import stats
from statsmodels.stats.multitest import multipletests

# %%
import seaborn as sns

# %%
# Define input and output files
if 'snakemake' in locals():
    tissue = snakemake.wildcards[0]

    sig_thresh = snakemake.params[0]
    cellprop_cutoff = snakemake.params[1]

    adata_fp = snakemake.input[0]
    inter_fp = snakemake.input[1]
    pathways_fp = snakemake.input[2]
    para_pathways_fp = snakemake.input[3]
    TFs_fp = snakemake.input[4]
    ctProp_fp = snakemake.input[5]

    output1_fp = snakemake.output[0]
    output2_fp = snakemake.output[1]

else:
    tissue = 'brain'
    
    sig_thresh = 0.05
    cellprop_cutoff = 0.05

    adata_fp = 'results/ST/{0}_wImages.h5ad'.format(tissue)
    inter_fp = 'results/Misty/{0}/CTpathways_diffInteractions.csv'.format(tissue)
    pathways_fp = 'results/ST/functional/{0}_activities_pathways.csv'.format(tissue)
    para_pathways_fp = 'results/ST/Misty/{0}/CTpathways/paraviews.csv'.format(tissue)
    TFs_fp = 'results/ST/functional/{0}_activities_GRNs.csv'.format(tissue)
    ctProp_fp = 'results/ST/ST_{0}_deconvoluted.csv'.format(tissue)

activities_fp ={'pathways': pathways_fp, 'para_pathways': para_pathways_fp, 'TFs': TFs_fp, 'cell_prop': ctProp_fp}

# %%
adata = sc.read_h5ad(adata_fp)
del adata.layers['SCT']
adata

# %%
interactions = pd.read_csv(inter_fp, index_col=0)
interactions = interactions[interactions['p.adj'] < sig_thresh].reset_index(drop=True)

# %%
for name, path in activities_fp.items():
    data = pd.read_csv(path, index_col=0).loc[adata.obs.index]
    data.columns = [re.sub("-", "", func) for func in data.columns]
    adata.obsm[name] = data

# %%
corrs = []
ttests = []
for index, row in interactions.iterrows():
    ad = adata[adata.obsm['cell_prop'][row['Target']] >= cellprop_cutoff].copy()
    print(ad.n_obs)

    # put together predictor and TF activities
    data = pd.concat([ad.obs.filter(['library_id', 'condition'], axis = 1),\
    ad.obsm[['para_pathways' if row['view'] == 'para' else 'pathways'][0]].filter([row['Predictor']], axis = 1),\
    ad.obsm['TFs']
        ], axis = 1)

    # computing the correlation
    corr = data.groupby(['library_id', 'condition']).corr().iloc[0::(ad.obsm['TFs'].shape[1] + 1), 1:]
    corr = corr.reset_index().rename({'level_2': 'Predictor'}, axis = 1)
    corr['inter'] = row['view'] + ':' + row['Predictor'] + ':' + row['Target']
    corrs.append(corr)

    temp = corr.drop(['library_id', 'Predictor', 'inter'], axis = 1)
    temp = pd.melt(temp, id_vars='condition', var_name='TF', value_name='Correlation')

    ttest = temp.groupby('TF').apply(lambda df: stats.ttest_ind(df['Correlation'][df['condition'] == 'flight'], df['Correlation'][df['condition'] == 'ground']))
    ttest = pd.DataFrame([[a, b] for a, b in ttest.values], columns = ['t_value', 'p_val'], index = ttest.index.values)
    _,ttest['padj'],_,_= multipletests(ttest['p_val'], alpha=sig_thresh, method = 'fdr_bh')
    ttest = ttest.reset_index().rename({'index': 'TF'}, axis = 1)
    ttest['Predictor'] = row['Predictor']
    ttest['inter'] = row['view'] + ':' + row['Predictor'] + ':' + row['Target']
    ttests.append(ttest)

ttests = pd.concat(ttests, axis = 0)
corrs = pd.concat(corrs, axis = 0)

#%%
if 'snakemake' in locals():
    ttests.to_csv(output1_fp, sep=',')
    corrs.to_csv(output2_fp, sep=',')

# # %%
# sig_associations = ttests[ttests['padj'] < sig_thresh].sort_values('padj', axis = 0, ascending = True).reset_index(drop=True)

# # %%
# mice = adata.obs.filter(['library_id', 'mouse']).drop_duplicates().reset_index(drop=True)

# # %%
# fig, axs = plt.subplots(sig_associations.shape[0], 5 , figsize=(25, sig_associations.shape[0] * 6))

# for index, row in sig_associations.iterrows():
#     print(row)
#     toplot = corrs[corrs['inter'] == row['inter']].filter(['library_id', 'condition', 'Predictor', row['TF']], axis = 1)
#     toplot = pd.merge(toplot, mice, how='left', on = ['library_id']).sort_values('mouse', axis = 0)
#     toplot['condition'] = toplot['condition'].str.capitalize()

#     flight_mouse = toplot[toplot['condition'] == 'Flight'].sort_values(row['TF'], axis = 0, ascending=[False if row['t_value'] > 0 else True][0])['library_id'].to_list()[0]
#     ground_mouse = toplot[toplot['condition'] == 'Ground'].sort_values(row['TF'], axis = 0, ascending=[True if row['t_value'] > 0 else False][0])['library_id'].to_list()[0]

#     sns.boxplot(x = "condition", y = row['TF'], data = toplot, ax = axs[index, 0], color = 'lightgrey', fliersize = 0, order = ['Flight', 'Ground'], width = 0.6)
#     sns.stripplot(x = "condition", y = row['TF'], data = toplot, hue = 'mouse', ax = axs[index, 0], order = ['Flight', 'Ground'])
#     axs[index, 0].set_title(ct_annot[ct_annot['clusterID'] == row['inter'].split(':')[-1]].iloc[0,2] + ': '+ row['inter'].split(':')[0] + row['Predictor'] + ' ~ ' + row['TF'])
#     axs[index, 0].set_ylabel('Correlation')

#     x1, x2 = 0, 1
#     y, h, col = toplot[row['TF']].max() * (1.1), toplot[row['TF']].max() * (0.1), 'k'
#     axs[index, 0].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

#     pval_text = '{:.2e}'.format(row['padj'])
#     if row['padj'] > sig_thresh: 
#         pval_text = 'ns = ' + pval_text
#     else: 
#         pval_text = 'sig = ' + pval_text
#     axs[index, 0].text((x1+x2)*.5, y+h, pval_text, ha='center', va='bottom', color=col)

#     for mouse, plot_column in zip([flight_mouse, ground_mouse], [1,3]):
#         temp = adata.copy()
#         temp = temp[temp.obs.library_id == mouse, :]
#         if row['inter'].split(':')[0] == 'para':
#             plotting = 'para_pathways'
#         else:
#             plotting = 'pathways'

#         set_na = temp.obsm['cell_prop'][row['inter'].split(':')[-1]] < cellprop_cutoff
#         temp.obsm[plotting][set_na] = np.NaN
#         temp.obsm['TFs'][set_na] = np.NaN

#         vmax = temp.obsm[plotting][row['Predictor']].quantile(0.98)
#         vmin = temp.obsm[plotting][row['Predictor']].quantile(0.02)

#         acts = dc.get_acts(temp, plotting)

#         sc.pl.spatial(acts, img_key=None, library_id=mouse, colorbar_loc=None,\
#             color=row['Predictor'], size=1.5, na_color = '#A69F9F', legend_loc=None, show=False, ax=axs[index, plot_column], vmax = vmax, vmin = vmin)
#         axs[index, plot_column].set_title(acts.obs['mouse'][0] + ': ' + row['inter'].split(':')[0] + row['Predictor'])
#         axs[index, plot_column].set_facecolor('#D9D9D9')
#         axs[index, plot_column].set_ylabel('')
#         axs[index, plot_column].set_xlabel('')

#         vmax = temp.obsm['TFs'][row['TF']].quantile(0.98)
#         vmin = temp.obsm['TFs'][row['TF']].quantile(0.02)

#         acts = dc.get_acts(temp, 'TFs')

#         sc.pl.spatial(acts, img_key=None, library_id=mouse, colorbar_loc=None,\
#             color=row['TF'], size=1.5, na_color = '#A69F9F', legend_loc=None,  show=False, ax=axs[index, plot_column + 1], vmax = vmax, vmin = vmin)
#         axs[index, plot_column + 1].set_title(acts.obs['mouse'][0] + ': ' + row['TF'])
#         axs[index, plot_column + 1].set_facecolor('#D9D9D9')
#         axs[index, plot_column + 1].set_ylabel('')
#         axs[index, plot_column + 1].set_xlabel('')

#         plt.tight_layout()

# # %%


# # %%


# # %%



