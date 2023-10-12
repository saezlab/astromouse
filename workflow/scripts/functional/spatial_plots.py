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

    adata_fp = snakemake.input[0]
    functional_fp = snakemake.input[1]

    output_fp = snakemake.output[0]

else:
    tissue = 'brain'

    adata_fp = 'results/ST/{0}_wImages.h5ad'.format(tissue)
    functional_fp = 'results/ST/functional/{0}_activities_pathways.csv'.format(tissue)



# %%
#Load anndata object
adata = sc.read_h5ad(adata_fp)
del adata.layers['SCT']
adata

# %%
#Load spot-level data to plot
activities = pd.read_csv(functional_fp, index_col=0, sep=',')
activities = activities.loc[adata.obs.index,:]
activities.columns = [re.sub("-", "", func) for func in activities.columns]

#Add functional data to anndata
adata.obsm['acts'] = activities
acts = dc.get_acts(adata, 'acts') #Reformat anndata
print(acts)

# %%
#Extract upper and lower data limits across samples for each feature
lims = pd.DataFrame({ 'llim' : [np.min(acts.X[:,ii]) for ii in range(acts.n_vars)], 'ulim': [np.max(acts.X[:,ii]) for ii in range(acts.n_vars)]}, index = acts.var_names.values)
lims['lim'] = [np.max(abs(acts.X[:,ii])) for ii in range(acts.n_vars)]
print('Max and min values per pathway')
print(lims)


if tissue == 'brain':
    with PdfPages(output_fp) as output_pdf:
        #Loop over features (i.e. called pathway)
        for pathway in acts.var.index.values:
            fig, axs = plt.subplots(3, 4, figsize=(23, 15))
            axs = axs.flatten()

            #Loop over samples
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
        #Loop over samples
        for library in acts.obs.filter(['library_id','mouse'], axis = 1).drop_duplicates().sort_values('library_id')['library_id']:
            #[os.path.basename(os.path.dirname(sample)) for sample in sample_paths]
            fig, axs = plt.subplots(4, 4, figsize=(23, 20))
            axs = axs.flatten()

            ad = acts[acts.obs.library_id == library, :]#.copy()

            #Loop over features (i.e. called pathway)
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


