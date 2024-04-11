# %%
import scanpy as sc
import pandas as pd
import decoupler as dc
from glob import glob
import re
import numpy as np
import os

# %%
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# %%
if 'snakemake' in locals():
    tissue = snakemake.wildcards[0]
    view_type = snakemake.wildcards[1]

    adata_fp = snakemake.input.get('data')

    pathways_fp = snakemake.input.get('pathways')
    TFs_fp = snakemake.input.get('GRNs')
    cellprop_fp = snakemake.input.get('cellprops')
    corr_fp = snakemake.input.get('correlations')
    paraviews_fp = snakemake.input.get('paraviews')
    cell_type_annot_fp = snakemake.input.get('cell_annot')

    
    importances_fp = [snakemake.input.get('importances')]
    interactions_fp = [snakemake.input.get('diffInteractions')]

    significance_threshold = snakemake.params[0]
    cellprop_cutoff = snakemake.params[1]

    output_files = snakemake.output

else:
    tissue = 'brain'
    view_type = 'celltype'

    adata_fp = 'results/ST/{0}_wImages.h5ad'.format(tissue)
    importances_fp = sorted(glob('results/Misty/{0}/*_importances.csv'.format(tissue)))
    interactions_fp = sorted(glob('results/Misty/{0}/*_diffInteractions.csv'.format(tissue)))
    pathways_fp = 'results/ST/functional/{0}_activities_pathways.csv'.format(tissue)
    TFs_fp = 'results/ST/functional/{0}_activities_TFs.csv'.format(tissue)
    cellprop_fp = 'results/ST/ST_{0}_deconvoluted.csv'.format(tissue)
    corr_fp = 'results/ST/Misty/{0}/{1}_Corr.csv'.format(tissue, view_type)
    cell_type_annot_fp ='data/original/MO/MO_cluster_metadata.csv'

    interactions_fp = interactions_fp[0:1]
    importances_fp = importances_fp[0:1]

    output_files = 'celltype.test'

    significance_threshold = 0.05
    cellprop_cutoff = 0.05

functional_fps = [pathways_fp, TFs_fp, cellprop_fp]
functional_fps = [f for f in functional_fps if f is not None]

correlations = pd.read_csv(corr_fp, index_col=0, sep=',')
correlations['Interaction'] = correlations[['view', 'Predictor', 'Target']].agg('_'.join, axis=1)
correlations = correlations.filter(['sample', 'corr', 'Interaction']).rename({'corr': 'Correlation'}, axis = 1)


# %%
models = [os.path.basename(fp).split('_')[0] for fp in interactions_fp]

# %%
adata = sc.read_h5ad(adata_fp)
del adata.layers['SCT']
adata

# %%

# dfs = []
# for file in functional_fps:
#     if 'deconvoluted' in file:
#         temp = pd.read_csv(file, index_col=0, sep=',').where(temp >= cellprop_cutoff)
#     else:
#         temp = pd.read_csv(file, index_col=0, sep=',')
    
#     dfs.append(temp)

dfs = [pd.read_csv(file, index_col=0, sep=',') for file in functional_fps]
activities = pd.concat(dfs, join='outer', axis=1)
activities = activities.loc[adata.obs.index,:]
activities.columns = [re.sub("-", "", func) for func in activities.columns]

adata.obsm['acts'] = activities
acts = dc.get_acts(adata, 'acts')
print(acts)

# %%
#Read in paraview and associate with corresponding cells
paraview = pd.read_csv(paraviews_fp, index_col=0, sep = ',').filter(adata.obs.index, axis = 0)
acts.obsm['paraview'] = paraview

# read in celltype cluster mapping
if cell_type_annot_fp is not None:
    cell_annot = pd.read_csv(cell_type_annot_fp, sep = ',')

# %%
df = [pd.read_csv(file, index_col=0, sep=',') for file in interactions_fp]
df = pd.concat(df, join='outer', axis=0)
mask = pd.DataFrame({'target':[ (target in activities.columns) for target in df['Target'] ],\
    'pred': [ (predictor in activities.columns) for predictor in df['Predictor'] ],\
    'sig': (df['p.adj'] <= significance_threshold)})
if mask['sig'].sum() == 0:
    df = df.sort_values('p.adj')[mask.filter(['target', 'pred']).all(axis='columns')].head(5).reset_index(drop=True)
else:
    df = df[mask.all(axis= 'columns')]
print('Significant interactions: ', df.shape[0])

# %%
interactions = {}
for model in models:
    temp = df[df['model']==model]
    if temp.shape[0] >0:
        interactions[model] = temp

# %%
df = [pd.read_csv(file, index_col=0, sep=',') for file in importances_fp]
df = pd.concat(df, join='outer', axis=0)

# %%
importances = {}
for model in models:
    temp = df[df['model']==model]
    if temp.shape[0] >0:
        importances[model] = temp

# %%
for key in interactions.keys():
    with PdfPages(output_files[0]) as output_pdf:
        current_inter = interactions[key].copy()
        current_inter['group'] = current_inter['view'] + '_' + current_inter['Predictor']
        current_inter = current_inter.sort_values('group', axis = 0)

        current_imp = importances[key].copy()
        current_imp = pd.merge(current_imp, correlations, how='left', on = ['sample', 'Interaction'])

        # select predictors to plot together
        ordered_predictors = current_inter.groupby(['view','group']).agg({'p.adj' : min}).sort_values(['view', 'p.adj']).reset_index().groupby('view').head(5)['group'].to_list()

        for current_pred in ordered_predictors:
            inter_to_plot = current_inter[current_inter['group'] == current_pred].head(5).copy().reset_index()
            inter_to_plot['inter'] = inter_to_plot['view'] + '_' + inter_to_plot['Predictor'] + '_' + inter_to_plot['Target']

            fig, axs = plt.subplots(inter_to_plot.shape[0], 6, figsize=(25, inter_to_plot.shape[0] * 5))

            mice = []

            #choosing which mice to plot
            #make boxplots
            for index, row in inter_to_plot.iterrows():
                imp_to_plot = current_imp[current_imp['Interaction'] == row['inter']].copy()
                imp_to_plot = imp_to_plot.sort_values('mouse')

                predictor_name = row['Predictor']
                target_name = row['Target']
                if view_type in ['celltype', 'CTpathways']:
                    target_name = cell_annot[cell_annot['clusterID'] == target_name]['clusterAbrv'].values[0]
                elif view_type in ['celltype', 'pathwaysCT']:
                    predictor_name = cell_annot[cell_annot['clusterID'] == predictor_name]['clusterAbrv'].values[0]

                flight_mouse = imp_to_plot[imp_to_plot['condition'] == 'Flight']
                ground_mouse = imp_to_plot[imp_to_plot['condition'] == 'Control']

                if np.sign(row['t.value']) > 0:
                    asc = False
                else:
                    asc = True

                flight_mouse = flight_mouse.sort_values('Importance', ascending = asc).head(1)
                ground_mouse = ground_mouse.sort_values('Importance', ascending = (not asc)).head(1)
                
                idx = [0 if inter_to_plot.shape[0] == 1 else (index, 0)]
                sns.boxplot(x = "condition", y = "Importance", data = imp_to_plot, ax = axs[idx[0]],\
                    color = 'lightgrey', fliersize = 0, order = ['Flight', 'Control'], width = 0.6)

                x1, x2 = 0, 1
                y, h, col = imp_to_plot["Importance"].max() * (1.1), imp_to_plot["Importance"].max() * (0.1), 'k'
                axs[idx[0]].plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)

                pval_text = '{:.2e}'.format(row['p.adj'])
                if row['p.adj'] > significance_threshold: 
                    pval_text = 'ns = ' + pval_text
                else: 
                    pval_text = 'p.adj = ' + pval_text
                axs[idx[0]].text((x1+x2)*.5, y+h, pval_text, ha='center', va='bottom', color=col)

                if tissue == 'brain':
                    sns.stripplot(x = "condition", y = "Importance", data = imp_to_plot, hue = 'mouse', ax = axs[idx[0]], order = ['Flight', 'Control'])
                else:
                    sns.stripplot(x = "condition", y = "Importance", data = imp_to_plot, ax = axs[idx[0]], order = ['Flight', 'Control'])
                axs[idx[0]].set_title(predictor_name + ' -> ' + target_name)

                #plot correlations
                idx = [5 if inter_to_plot.shape[0] == 1 else (index, 5)]
                sns.boxplot(x = "condition", y = "Correlation", data = imp_to_plot, ax = axs[idx[0]],\
                    color = 'lightgrey', fliersize = 0, order = ['Flight', 'Control'], width = 0.6)

                if tissue == 'brain':
                    sns.stripplot(x = "condition", y = "Correlation", data = imp_to_plot, hue = 'mouse', ax = axs[idx[0]], order = ['Flight', 'Control'])
                else:
                    sns.stripplot(x = "condition", y = "Correlation", data = imp_to_plot, ax = axs[idx[0]], order = ['Flight', 'Control'])
                axs[idx[0]].set_title(predictor_name + ' -> ' + target_name)

                mice.append(flight_mouse)
                mice.append(ground_mouse)

            #pool mice to plot
            meta = pd.merge(pd.concat(mice, join='outer', axis=0), inter_to_plot, how = 'left', left_on='Interaction', right_on='inter')
            func_to_plot = list(set(meta['Predictor'].to_list() + meta['Target'].to_list())) #which features will be plotted

            #identify limits for colorbars
            lims = pd.DataFrame({ 'llim' : [np.min(acts[:,f].X) for f in func_to_plot],\
            'ulim': [np.max(acts[:,f].X) for f in func_to_plot]},\
                index = func_to_plot)
            lims['lim'] = [np.max(abs(acts[:,f].X)) for f in func_to_plot]

            # iterate over rows (interactions) of the plot
            for index, row in inter_to_plot.iterrows():
                mice_to_plot = meta[meta['Interaction']  == row['inter']].sort_values('condition', ascending=False)
                
                for (_, mouse), plot_column in zip(mice_to_plot.iterrows(), [1,3]):

                    temp = acts.copy()
                    if inter_to_plot['view'].to_list()[0] == 'para':
                        plotting = 'paraview'
                    else:
                        plotting = 'acts'
                        # masking cells below a specific threshold (only in intra/intra_act)
                        if (view_type == 'celltype') or (view_type == 'pathwaysCT'):
                            temp.obsm['acts'] = temp.obsm['acts'].where(temp.obsm['acts'] >= cellprop_cutoff)

                    temp = dc.get_acts(temp, plotting)
                    temp = temp[temp.obs.library_id == mouse['sample'], :] #select cells of one mouse

                    idx = [plot_column if inter_to_plot.shape[0] == 1 else (index, plot_column)]
                    sc.pl.spatial(temp, img_key=None, library_id=mouse['sample'],\
                        color=mouse['Predictor'], size=1.5, na_color = '#A69F9F', legend_loc=None, show=False, ax=axs[idx[0]]) #vmin = (lims.loc[mouse['Predictor'], 'llim']*1.1),\vmax = (lims.loc[mouse['Predictor'], 'ulim']*1.1)
                    axs[idx[0]].set_title(mouse['mouse'] + ': ' + predictor_name)
                    axs[idx[0]].set_facecolor('#D9D9D9')
                    axs[idx[0]].set_ylabel('')
                    axs[idx[0]].set_xlabel('')

                    temp = acts.copy()
                    # masking cells below a specific threshold (only in intra/intra_act)
                    if (view_type == 'celltype') or (view_type == 'CTpathways'):
                        temp.obsm['acts'] = temp.obsm['acts'].where(temp.obsm['acts'] >= cellprop_cutoff)
                        temp = dc.get_acts(temp, 'acts')
                    temp = temp[temp.obs.library_id == mouse['sample'], :]

                    idx = [plot_column  + 1 if inter_to_plot.shape[0] == 1 else (index, plot_column + 1)]
                    sc.pl.spatial(temp, img_key=None, library_id=mouse['sample'],\
                        color=mouse['Target'], size=1.5, na_color = '#A69F9F', legend_loc=None, show=False, ax=axs[idx[0]]) #vmin = (lims.loc[mouse['Target'], 'llim']*1.1),\vmax = (lims.loc[mouse['Target'], 'ulim']*1.1), 
                    axs[idx[0]].set_title(mouse['mouse'] + ': ' + target_name)
                    axs[idx[0]].set_facecolor('#D9D9D9')
                    axs[idx[0]].set_ylabel('')
                    axs[idx[0]].set_xlabel('')
                
            plt.suptitle('Diff. interactions in ' + inter_to_plot['view'].to_list()[0] + ' for ' + current_pred.split('_')[-1], fontsize = 15)
            plt.tight_layout()
            output_pdf.savefig(fig)
            plt.close()
                
