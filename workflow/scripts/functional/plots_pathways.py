import scanpy as sc
import pandas as pd
import decoupler as dc
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Define input and output files
if 'snakemake' in locals():
    adata_fp = snakemake.input[0]
    output_fp = snakemake.output[0]
    tissue = snakemake.wildcards[0]
    normalisation = snakemake.params[0]
    top_genes = snakemake.params[1]


adata_spatial = sc.read_h5ad(adata_fp)

if normalisation == 'log1p':
    sc.pp.normalize_total(adata_spatial, inplace=True)
    sc.pp.log1p(adata_spatial)
elif normalisation == 'SCT':
    adata_spatial.layers['counts'] = adata_spatial.X
    adata_spatial.X = adata_spatial.layers['SCT'].copy()
    del adata_spatial.layers['SCT']
    adata_spatial


# fig, axs = plt.subplots(3, 4, figsize=(23, 15))
# axs = axs.flatten()

# for i, library in enumerate(
#     adata_spatial.obs.filter(['library_id','mouse'], axis = 1).drop_duplicates().sort_values('mouse')['library_id']
#     #[os.path.basename(os.path.dirname(sample)) for sample in sample_paths]
# ):
#     ad = adata_spatial[adata_spatial.obs.library_id == library, :]#.copy()
#     sc.pl.spatial(
#         ad,
#         img_key=None,#"hires",
#         library_id=library,
#         color="seurat_clusters",
#         size=1.5,
#         palette=sc.pl.palettes.default_20,
#         legend_loc=None,
#         show=False,
#         ax=axs[i],
#     )

# plt.tight_layout()

model = dc.get_progeny(organism='mouse', top=top_genes)
dc.run_mlm(mat=adata_spatial, net=model, source='source', target='target', weight='weight', verbose=True, use_raw=False)

acts = dc.get_acts(adata_spatial, obsm_key='mlm_estimate')
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
                    vmin = (lims.loc[pathway, 'llim']*1.1),
                    vmax = (lims.loc[pathway, 'ulim']*1.1),
                    color_map = 'BrBG',
                    vcenter = 0,
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
                    vmin = (lims.loc[pathway, 'llim']*1.1),
                    vmax = (lims.loc[pathway, 'ulim']*1.1),
                    color_map = 'BrBG',
                    vcenter = 0,
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