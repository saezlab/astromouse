import scanpy as sc
import pandas as pd
import decoupler as dc
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Define input and output files
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
        model = pd.read_csv(os.path.join(net_fp, tissue + '.csv'), sep=',', index_col=0)
else:
    raise ValueError('The "network" wildcard can only take on "pathways", "GRNS" or "TFs" as value, to run either Progeny, celloracle GRNs or Dorothea regulons')


# %%
dc.decouple(mat=adata, net=model, source='source', target='target', weight='weight', methods = conf.get('method'),  verbose=True, use_raw=False)

# %%
print(adata)

acts = dc.get_acts(adata, obsm_key='{0}_estimate'.format(conf.get('method')))
print(acts)

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
                    # vmin = 0,
                    # vmax = 17589,
                    # vmin = (lims.loc[pathway, 'llim']*1.1),
                    # vmax = (lims.loc[pathway, 'ulim']*1.1),
                    color_map = 'BrBG',
                    colorbar_loc = None,
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
                    color_map = 'BrBG',
                    colorbar_loc = None,
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