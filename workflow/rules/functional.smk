#input function for activity prediction
def activities_inputs(wildcards):
    files = {'data': 'results/ST/{wildcards.tissue}_wImages.h5ad'.format(wildcards=wildcards)}
    if wildcards.network == 'GRNs':
        files['net'] = 'results/MO/celloracle/{wildcards.tissue}/GRNs/'.format(wildcards=wildcards)
    return files

#infer activity for Visium data
#based on GRN (for TFs) or the progeny resource (for pathways)
rule get_PDactivities:
    input:
        unpack(activities_inputs)
    output:
        act = 'results/ST/functional/{tissue}_activities_{network}.csv' #network being either 'pathways' or 'TFs'
    params:
        lambda w: config['functional'][w.network]
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/functional/compute_activities.py"

#spatial plot of pathway activities
rule plot_pathways:
    input:
        adata = 'results/ST/{tissue}_wImages.h5ad',
        functional = 'results/ST/functional/{tissue}_activities_pathways.csv'
    output:
        'plots/functional/{tissue}/pathways.pdf'
    conda:
        "../envs/astromouse.yml"
    script:
        '../scripts/functional/spatial_plots.py'

#spatial plot of pathway paraviews (feature in Misty models)
rule plot_pathway_paraviews:
    input:
        adata = 'results/ST/{tissue}_wImages.h5ad',
        functional = 'results/ST/Misty/{tissue}/CTpathways/paraviews.csv'
    output:
        'plots/functional/{tissue}/pathway_paraviews.pdf'
    conda:
        "../envs/astromouse.yml"
    script:
        '../scripts/functional/spatial_plots.py'

#spatial plot of cell type abundances (from stereoscope)
rule plot_spatial_stereoscope:
    input:
        adata = 'results/ST/{tissue}_wImages.h5ad',
        functional = 'results/ST/ST_{tissue}_deconvoluted.csv',
        annotations = 'data/MO_cluster_metadata.csv'
    params:
        cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff')
    output:
        'plots/functional/{tissue}/stereoscope.pdf'
    conda:
        "../envs/astromouse.yml"
    script:
        '../scripts/functional/stereoscope_plots.py'

#plot showing abundance of celltypes per sample
rule plot_stereoscope:
    input:
        adata = 'results/ST/{tissue}_wImages.h5ad',
        functional = 'results/ST/ST_{tissue}_deconvoluted.csv',
        annotations = 'data/MO_cluster_metadata.csv'
    params:
        cellprop_cutoff = 0.1
    output:
        'plots/functional/{tissue}/stereoscope_proportions.pdf'
    conda:
        "../envs/astromouse.yml"
    script:
        '../scripts/functional/stereoscope_proportions.py'

#spatial plot of sample images
rule plot_slides:
    input:
        adata = 'results/ST/{tissue}_wImages.h5ad'
    output:
        'plots/functional/{tissue}/slides.pdf'
    conda:
        "../envs/astromouse.yml"
    script:
        '../scripts/functional/slide_plots.py'

#not used in analysis
# rule plot_ROIs:
#     input:
#         adata = 'results/ST/{tissue}_wImages.h5ad'
#     output:
#         'plots/functional/{tissue}/spatial_clusters.pdf',
#         'plots/functional/{tissue}/ROIs.pdf'
#     conda:
#         "../envs/astromouse.yml"
#     script:
#         '../scripts/functional/plot_ROIs.py'