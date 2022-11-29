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
# Define input and output files
if 'snakemake' in locals():
    tissue = snakemake.wildcards[0]
    cellprop_cutoff = snakemake.params[0]

    adata_fp = snakemake.input[0]
    functional_fp = snakemake.input[1]
    annot_fp = snakemake.input[2]

    output_fp = snakemake.output[0]

else:
    tissue = 'brain'
    cellprop_cutoff = 0.05


    adata_fp = 'results/ST/{0}_wImages.h5ad'.format(tissue)
    functional_fp = 'results/ST/ST_{0}_deconvoluted.csv'.format(tissue)
    annot_fp = 'data/original/MO/MO_cluster_metadata.csv'

# %%
annot = pd.read_csv(annot_fp)


# %%
adata = sc.read_h5ad(adata_fp)
del adata.layers['SCT']
adata

# %%
activities = pd.read_csv(functional_fp, index_col=0, sep=',')
activities = activities.loc[adata.obs.index,:]
activities.columns = annot['clusterAbrv']
activities = activities.where(activities >= cellprop_cutoff)

adata.obsm['acts'] = activities.reindex(sorted(activities.columns), axis=1)
acts = dc.get_acts(adata, 'acts')
print(acts)

# %%
lims = pd.DataFrame({ 'llim' : [np.min(acts.X[:,ii]) for ii in range(acts.n_vars)], 'ulim': [np.max(acts.X[:,ii]) for ii in range(acts.n_vars)]}, index = acts.var_names.values)
lims['lim'] = [np.max(abs(acts.X[:,ii])) for ii in range(acts.n_vars)]
print('Max and min values per pathway')
print(lims)


if tissue == 'brain':
    with PdfPages(output_fp) as output_pdf:
        for pathway in acts.var.index.values:
            fig, axs = plt.subplots(3, 4, figsize=(23, 15))
            axs = axs.flatten()

            for i, library in enumerate(
                acts.obs.filter(['library_id','mouse'], axis = 1).drop_duplicates().sort_values('mouse')['library_id']
                #[os.path.basename(os.path.dirname(sample)) for sample in sample_paths]
            ):
                ad = acts[acts.obs.library_id == library, :]#.copy()
                sc.pl.spatial(
                    ad,
                    img_key=None,
                    library_id=library,
                    color=pathway,
                    size=1.5,
                    legend_loc=None,
                    show=False,
                    na_color = '#A69F9F',
                    # vmin = 0,
                    # vmax = 17589,
                    # vmin = (lims.loc[pathway, 'llim']*1.1),
                    # vmax = (lims.loc[pathway, 'ulim']*1.1),
                    # color_map = 'BrBG',
                    # colorbar_loc = None,
                    # vcenter = 0,
                    ax=axs[i],
                )
                axs[i].set_title(ad.obs['mouse'][0])
                axs[i].set_facecolor('#D9D9D9')
                axs[i].set_ylabel('')
                axs[i].set_xlabel('')

            plt.suptitle(pathway, fontsize = 15)
            plt.tight_layout()
            output_pdf.savefig(fig)
            plt.close()

if tissue == 'heart':
    with PdfPages(output_fp) as output_pdf:
        for library in acts.obs.filter(['library_id','mouse'], axis = 1).drop_duplicates().sort_values('library_id')['library_id']:
            #[os.path.basename(os.path.dirname(sample)) for sample in sample_paths]
            fig, axs = plt.subplots(4, 4, figsize=(23, 20))
            axs = axs.flatten()

            ad = acts[acts.obs.library_id == library, :]#.copy()

            for i, pathway in enumerate(acts.var.index.values):

                sc.pl.spatial(
                    ad,
                    img_key=None,
                    library_id=library,
                    color=pathway,
                    size=1.5,
                    legend_loc=None,
                    show=False,
                    # vmin = 0,
                    # vmax = 17589,
                    # vmin = (lims.loc[pathway, 'llim']*1.1),
                    # vmax = (lims.loc[pathway, 'ulim']*1.1),
                    # color_map = 'BrBG',
                    # colorbar_loc = None,
                    # vcenter = 0,
                    ax=axs[i],
                )
                axs[i].set_title(pathway)
                axs[i].set_facecolor('#D9D9D9')
                axs[i].set_ylabel('')
                axs[i].set_xlabel('')

            for ax in axs[lims.shape[0]:]:
                ax.set_axis_off()

            plt.suptitle(ad.obs['mouse'][0], fontsize = 15)
            plt.tight_layout()
            output_pdf.savefig(fig)
            plt.close()


