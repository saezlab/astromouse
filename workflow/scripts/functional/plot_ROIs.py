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

    output1_fp = snakemake.output[0]
    output2_fp = snakemake.output[1]

else:
    tissue = 'brain'

    adata_fp = 'results/ST/{0}_wImages.h5ad'.format(tissue)



# %%
adata = sc.read_h5ad(adata_fp)
del adata.layers['SCT']
adata

# %%
adata.obs.seurat_clusters = adata.obs.seurat_clusters.astype('int')

# %%
with PdfPages(output1_fp) as output_pdf:
    fig, axs = plt.subplots(3, 4, figsize=(23, 15))
    axs = axs.flatten()

    for i, library in enumerate(adata.obs.filter(['library_id','mouse'], axis = 1).drop_duplicates().sort_values('mouse')['library_id']):
        ad = adata[(adata.obs.library_id == library), :]#.copy()
        sc.pl.spatial(
            ad,
            # img_key=None,
            library_id=library,
            color='annot',
            # size=1.5,
            legend_loc=None,
            show=False,
            # vmin = 0,
            # vmax = 17589,
            # vmin = (lims.loc[pathway, 'llim']*1.1),
            # vmax = (lims.loc[pathway, 'ulim']*1.1),
            # color_map = 'BrBG',
            # colorbar_loc = None,
            # vcenter = 0,
            ax=axs[i]
        )
        axs[i].set_title(ad.obs['mouse'][0]+ ' ' + ad.obs['library_id'][0])
        # axs[i].set_facecolor('#3D3D3D')
        axs[i].set_ylabel('')
        axs[i].set_xlabel('')

    plt.tight_layout()
    output_pdf.savefig(fig)
    plt.close()


# %%
cluster_id = [[6, 8, 10, 11],[1, 9],[14],[16]]
cluster_name = ['Hippocampal neurons and dentate gyrus', 'Cortical layers', 'Caudate/putamen neurons', 'Choroid plexus and subventricular structures']

# %%
with PdfPages(output2_fp) as output_pdf:
    for clusters, name in zip(cluster_id, cluster_name):
        fig, axs = plt.subplots(3, 4, figsize=(23, 15))
        axs = axs.flatten()

        for i, library in enumerate(adata.obs.filter(['library_id','mouse'], axis = 1).drop_duplicates().sort_values('mouse')['library_id']):
            ad = adata[(adata.obs.library_id == library), :]#.copy()
            ad.obs.loc[[cl not in clusters for cl in ad.obs['seurat_clusters']], 'annot'] = np.NaN
            sc.pl.spatial(
                ad,
                # img_key=None,
                library_id=library,
                color='annot',
                # size=1.5,
                legend_loc=None,
                show=False,
                # vmin = 0,
                # vmax = 17589,
                # vmin = (lims.loc[pathway, 'llim']*1.1),
                # vmax = (lims.loc[pathway, 'ulim']*1.1),
                # color_map = 'BrBG',
                # colorbar_loc = None,
                # vcenter = 0,
                na_color = '#A69F9F',
                # palette = ['#F9D505'],
                ax=axs[i]
            )
            axs[i].set_title(ad.obs['mouse'][0]+ ' ' + ad.obs['library_id'][0])
            # axs[i].set_facecolor('#D9D9D9')
            axs[i].set_ylabel('')
            axs[i].set_xlabel('')

        
        plt.suptitle('ROI: ' + name)
        plt.tight_layout()
        output_pdf.savefig(fig)
        plt.close()
