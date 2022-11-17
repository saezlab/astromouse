############################
#   Set up and run misty models
#   sample by sample
############################

rule get_coords:
    input:
        data = 'results/ST/{tissue}_wImages.h5ad'
    output:
        coords = 'results/ST/Misty/{tissue}_coordinates.csv'
    conda:
        "../envs/astromouse.yml"
    resources:
        mem_mb=25000
    script:
        "../scripts/misty/extract_coordinates.py"

rule get_func_views:
    input:
        'results/ST/Misty/{tissue}_coordinates.csv',
        'results/ST/functional/{tissue}_activities_GRNs.csv',
        'results/ST/functional/{tissue}_activities_pathways.csv'
    params:
        skip = 'intra'
    output:
        view = 'results/ST/Misty/{tissue}/functional/views/{sample}_view.rds',
        paraview = 'results/ST/Misty/{tissue}/functional/views/{sample}_paraview.csv'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_celltype_views:
    input:
        'results/ST/Misty/{tissue}_coordinates.csv',
        'results/ST/ST_{tissue}_deconvoluted.csv'
    output:
        view = 'results/ST/Misty/{tissue}/celltype/views/{sample}_view.rds',
        paraview = 'results/ST/Misty/{tissue}/celltype/views/{sample}_paraview.csv'
    params:
        cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff')
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_pathwaysCT_views:
    input:
        'results/ST/Misty/{tissue}_coordinates.csv',
        'results/ST/functional/{tissue}_activities_pathways.csv',
        'results/ST/ST_{tissue}_deconvoluted.csv'
    output:
        view = 'results/ST/Misty/{tissue}/pathwaysCT/views/{sample}_view.rds',
        paraview = 'results/ST/Misty/{tissue}/pathwaysCT/views/{sample}_paraview.csv'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_CTpathways_views:
    input:
        'results/ST/Misty/{tissue}_coordinates.csv',
        'results/ST/ST_{tissue}_deconvoluted.csv',
        'results/ST/functional/{tissue}_activities_pathways.csv',
    params:
        cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff')
    output:
        view = 'results/ST/Misty/{tissue}/CTpathways/views/{sample}_view.rds',
        paraview = temp('results/ST/Misty/{tissue}/CTpathways/views/{sample}_paraview.csv')
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_combine_paraviews:
    input:
        lambda w: expand('results/ST/Misty/{{tissue}}/{{view_type}}/views/{sample}_paraview.csv', sample = config['samples'][w.tissue])
    output:
        'results/ST/Misty/{tissue}/{view_type}/paraviews.csv'
    shell:
        "awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} >> {output}"

rule run_views:
    input:
        view = 'results/ST/Misty/{tissue}/{view_type}/views/{sample}_view.rds'
    output: 
        directory('results/ST/Misty/{tissue}/{view_type}/models/{sample}')
    params:
        seed = config['misty'].get("random_seed", 42),
        bypass_intra = lambda wildcards: config['misty'][wildcards.view_type].get('bypass_intra', False)
    conda:
        "../envs/misty.yml"
    threads: 6
    resources:
        mem_mb=25000,
        disk_mb=1000,
        time='12:00:00'
    script:
        "../scripts/misty/run.R"

rule plot_misty_results:
    input:
        'data/original/ST/metadata_visium_{tissue}.csv',
        'data/original/MO/MO_cluster_metadata.csv',
        lambda w: expand('results/ST/Misty/{{tissue}}/{{view_type}}/models/{sample}', sample = config['samples'][w.tissue])
    output: 
        'plots/Misty/{tissue}/{view_type}_misty.pdf',
        'plots/Misty/{tissue}/{view_type}_misty_Flight.pdf',
        'plots/Misty/{tissue}/{view_type}_misty_Control.pdf'

    params:
        lambda w: config['misty'][w.view_type]['plots'][w.tissue]
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/plot_model_results.R"

# ############################
# #   Brain region specific models
# ############################

# rule get_region_celltype_views:
#     input:
#         'results/ST/Misty/{tissue}_coordinates.csv',
#         'results/ST/ST_{tissue}_deconvoluted.csv'
#     output:
#         view = 'results/ST/Misty/{tissue}/{region}/celltype/views/{sample}_view.rds',
#         paraview = 'results/ST/Misty/{tissue}/{region}/celltype/views/{sample}_paraview.csv'
#     params:
#         cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff'),
#         clusters = lambda w: config['regions']['brain'][w.region]
#     conda:
#         "../envs/misty.yml"
#     script:
#         "../scripts/misty/make_region_views.R"

# rule run_region_views:
#     input:
#         view = 'results/ST/Misty/{tissue}/{region}/celltype/views/{sample}_view.rds'
#     output: 
#         directory('results/ST/Misty/{tissue}/{region}/celltype/models/{sample}')
#     params:
#         seed = config['misty'].get("random_seed", 42),
#         bypass_intra = lambda wildcards: config['misty'][wildcards.view_type].get('bypass_intra', False)
#     conda:
#         "../envs/misty.yml"
#     threads: 6
#     resources:
#         mem_mb=25000,
#         disk_mb=1000,
#         time='12:00:00'
#     script:
#         "../scripts/misty/run.R"

############################
#   Get interactions from misty models
#   that are only found in one of the conditions
############################

rule get_dif_interactions:
    input:
        'data/original/ST/metadata_visium_{tissue}.csv',
        lambda w: expand('results/ST/Misty/{{tissue}}/{{view_type}}/models/{sample}', sample = config['samples'][w.tissue])
    output: 
        'results/Misty/{tissue}/{view_type}_importances.csv',
        'results/Misty/{tissue}/{view_type}_diffInteractions.csv'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/get_differential_interactions.R"

rule get_interaction_corr:
    input:
        interactions = 'results/Misty/{tissue}/{view_type}_diffInteractions.csv',
        view = 'results/ST/Misty/{tissue}/{view_type}/views/{sample}_view.rds'
    output:
        corr = temp('results/ST/Misty/{tissue}/{view_type}/correlations/{sample}_Corr.csv')
    params:
        corr = 'pearson'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/get_interactions_corr.R"

rule combine_interaction_corr:
    input:
        lambda w: expand('results/ST/Misty/{{tissue}}/{{view_type}}/correlations/{sample}_Corr.csv', sample = config['samples'][w.tissue])
    output:
        corr = 'results/Misty/{tissue}/{view_type}_Corr.csv'
    shell:
        "awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} >> {output}"


def dif_interactions_inputs(wildcards):
    files = {'data': 'results/ST/{wildcards.tissue}_wImages.h5ad'.format(wildcards=wildcards),\
            'importances': 'results/Misty/{wildcards.tissue}/{wildcards.view_type}_importances.csv'.format(wildcards=wildcards),\
            'diffInteractions': 'results/Misty/{wildcards.tissue}/{wildcards.view_type}_diffInteractions.csv'.format(wildcards=wildcards),
            'correlations': 'results/Misty/{wildcards.tissue}/{wildcards.view_type}_Corr.csv'.format(wildcards=wildcards),\
            'paraviews': 'results/ST/Misty/{wildcards.tissue}/{wildcards.view_type}/paraviews.csv'.format(wildcards=wildcards)
            }
    if (wildcards.view_type in ['celltype', 'pathwaysCT', 'CTpathways']):
        files['cellprops'] = 'results/ST/ST_{wildcards.tissue}_deconvoluted.csv'.format(wildcards=wildcards)
        files['cell_annot'] = 'data/original/MO/MO_cluster_metadata.csv'

    if (wildcards.view_type in ['functional', 'pathwaysCT', 'CTpathways']):
        files['pathways'] = 'results/ST/functional/{wildcards.tissue}_activities_pathways.csv'.format(wildcards=wildcards)

    if (wildcards.view_type == 'functional'):
        files['GRNs'] = 'results/ST/functional/{wildcards.tissue}_activities_GRNs.csv'.format(wildcards=wildcards)

    return files

rule plot_dif_interactions:
    input:
        unpack(dif_interactions_inputs)
    params:
        sign= 0.05,
        cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff')
    output: 
        'plots/Misty/{tissue}/{view_type}_diffplots.pdf'
    conda:
        "../envs/astromouse.yml"
    resources:
        mem_mb=25000
    script:
        "../scripts/misty/plot_differential_interactions.py"

############################
#   Relating pathway activities to TF activity changes
#   Only for specific cell types
############################

rule TF_pathway_corr:
    input:
        data = 'results/ST/{tissue}_wImages.h5ad',
        diff_inter = 'results/Misty/{tissue}/CTpathways_diffInteractions.csv',
        pathways = 'results/ST/functional/{tissue}_activities_pathways.csv',
        pathway_paras = 'results/ST/Misty/{tissue}/CTpathways/paraviews.csv',
        TFs = 'results/ST/functional/{tissue}_activities_GRNs.csv',
        cellProp = 'results/ST/ST_{tissue}_deconvoluted.csv'
    params:
        sign= 0.05,
        cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff')
    output:
        'results/Misty/{tissue}/interactions_TFPathway_ttests.csv',
        'results/Misty/{tissue}/interactions_TFPathway_corrs.csv'
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/misty/pathway_TF_corrs.py"

rule TF_pathway_spatial_plots:
    input:
        data = 'results/ST/{tissue}_wImages.h5ad',
        diff_inter = 'results/Misty/{tissue}/CTpathways_diffInteractions.csv',
        pathways = 'results/ST/functional/{tissue}_activities_pathways.csv',
        pathway_paras = 'results/ST/Misty/{tissue}/CTpathways/paraviews.csv',
        TFs = 'results/ST/functional/{tissue}_activities_GRNs.csv',
        cellProp = 'results/ST/ST_{tissue}_deconvoluted.csv',
        ct_annot = 'data/original/MO/MO_cluster_metadata.csv',
        ttests = 'results/Misty/{tissue}/interactions_TFPathway_ttests.csv',
        corrs = 'results/Misty/{tissue}/interactions_TFPathway_corrs.csv'
    params:
        sign= 0.05,
        cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff')
    output:
        'plots/Misty/{tissue}/interactions_TFPathway.pdf'
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/misty/plot_TFPathway_corrs.py"