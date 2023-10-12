############################
#   Set up and run misty models
#   sample by sample
############################

# extract coordinates and metadata of each visium spot
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

#make Misty views with TF activities (from GRN, inferred with decoupler) as targets
#with pathway activities as predictors
#also exports the paraview values
# (not used in analysis)
# rule get_func_views:
#     input:
#         'results/ST/Misty/{tissue}_coordinates.csv',
#         'results/ST/functional/{tissue}_activities_GRNs.csv',
#         'results/ST/functional/{tissue}_activities_pathways.csv'
#     params:
#         skip = 'intra'
#     output:
#         view = 'results/ST/Misty/{tissue}/functional/views/{sample}_view.rds',
#         paraview = 'results/ST/Misty/{tissue}/functional/views/{sample}_paraview.csv'
#     conda:
#         "../envs/misty.yml"
#     script:
#         "../scripts/misty/make_views.R"

#make Misty views with cell type abundances (from Stereoscope) as targets
#with cell type abundances as predictors
#also exports the paraview values
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

#make Misty views with pathway activities (from progeny, inferred with decoupler) as targets
#with cell type abundances as predictors
#also exports the paraview values
# (not used in analysis)
# rule get_pathwaysCT_views:
#     input:
#         'results/ST/Misty/{tissue}_coordinates.csv',
#         'results/ST/functional/{tissue}_activities_pathways.csv',
#         'results/ST/ST_{tissue}_deconvoluted.csv'
#     output:
#         view = 'results/ST/Misty/{tissue}/pathwaysCT/views/{sample}_view.rds',
#         paraview = 'results/ST/Misty/{tissue}/pathwaysCT/views/{sample}_paraview.csv'
#     conda:
#         "../envs/misty.yml"
#     script:
#         "../scripts/misty/make_views.R"

#make Misty views with cell type abundances as targets
#with pathway activities (from progeny, inferred with decoupler) as predictors
#also exports the paraview values
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

#combine the paraview values of all samples to one csv file
rule get_combine_paraviews:
    input:
        lambda w: expand('results/ST/Misty/{{tissue}}/{{view_type}}/views/{sample}_paraview.csv', sample = config['samples'][w.tissue])
    output:
        'results/ST/Misty/{tissue}/{view_type}/paraviews.csv'
    shell:
        "awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} >> {output}"

#generic rule to run a Misty model on a (sample's) Misty view
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

#analytical plots of the Misty models pooled across samples and conditions
rule plot_misty_results:
    input:
        'data/metadata_visium_{tissue}.csv',
        'data/MO_cluster_metadata.csv',
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

############################
#   Get interactions from misty models
#   that are only found in one of the conditions
############################

#extract interactions that are condition specific
rule get_dif_interactions:
    input:
        'data/metadata_visium_{tissue}.csv',
        lambda w: expand('results/ST/Misty/{{tissue}}/{{view_type}}/models/{sample}', sample = config['samples'][w.tissue])
    output: 
        'results/Misty/{tissue}/{view_type}_importances.csv',
        'results/Misty/{tissue}/{view_type}_diffInteractions.csv'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/get_differential_interactions.R"

#from differential interactions, extract the correlation between target-predictor
#sample by sample
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

#combine the interaction correlation of all samples to one csv file
rule combine_interaction_corr:
    input:
        lambda w: expand('results/ST/Misty/{{tissue}}/{{view_type}}/correlations/{sample}_Corr.csv', sample = config['samples'][w.tissue])
    output:
        corr = 'results/Misty/{tissue}/{view_type}_Corr.csv'
    shell:
        "awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} >> {output}"

#input function for rule below
#assigns the right files for plotting the predictors
def dif_interactions_inputs(wildcards):
    files = {'data': 'results/ST/{wildcards.tissue}_wImages.h5ad'.format(wildcards=wildcards),\
            'importances': 'results/Misty/{wildcards.tissue}/{wildcards.view_type}_importances.csv'.format(wildcards=wildcards),\
            'diffInteractions': 'results/Misty/{wildcards.tissue}/{wildcards.view_type}_diffInteractions.csv'.format(wildcards=wildcards),
            'correlations': 'results/Misty/{wildcards.tissue}/{wildcards.view_type}_Corr.csv'.format(wildcards=wildcards),\
            'paraviews': 'results/ST/Misty/{wildcards.tissue}/{wildcards.view_type}/paraviews.csv'.format(wildcards=wildcards)
            }
    if (wildcards.view_type in ['celltype', 'pathwaysCT', 'CTpathways']):
        files['cellprops'] = 'results/ST/ST_{wildcards.tissue}_deconvoluted.csv'.format(wildcards=wildcards)
        files['cell_annot'] = 'data/MO_cluster_metadata.csv'

    if (wildcards.view_type in ['functional', 'pathwaysCT', 'CTpathways']):
        files['pathways'] = 'results/ST/functional/{wildcards.tissue}_activities_pathways.csv'.format(wildcards=wildcards)

    if (wildcards.view_type == 'functional'):
        files['GRNs'] = 'results/ST/functional/{wildcards.tissue}_activities_GRNs.csv'.format(wildcards=wildcards)

    return files

#for each differential interaction make
#boxplots of predictor importances
#spatial plots of target/predictor data for two samples
#boxplots of the interaction correlation
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

#extract correlation of pathways and TF activities for differentially active pathways
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

#for each differential interaction make
#spatial plots of pathway/TF activity for two samples
#boxplots of the pathway-TF correlation
rule TF_pathway_spatial_plots:
    input:
        data = 'results/ST/{tissue}_wImages.h5ad',
        diff_inter = 'results/Misty/{tissue}/CTpathways_diffInteractions.csv',
        pathways = 'results/ST/functional/{tissue}_activities_pathways.csv',
        pathway_paras = 'results/ST/Misty/{tissue}/CTpathways/paraviews.csv',
        TFs = 'results/ST/functional/{tissue}_activities_GRNs.csv',
        cellProp = 'results/ST/ST_{tissue}_deconvoluted.csv',
        ct_annot = 'data/MO_cluster_metadata.csv',
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