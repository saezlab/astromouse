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
        view = 'results/ST/Misty/{tissue}/{sample}/functional_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_celltype_views:
    input:
        'results/ST/Misty/{tissue}_coordinates.csv',
        'results/ST/ST_{tissue}_deconvoluted.csv'
    output:
        view = 'results/ST/Misty/{tissue}/{sample}/celltype_view.rds'
    params:
        cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff')
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_patwayCT_views:
    input:
        'results/ST/Misty/{tissue}_coordinates.csv',
        'results/ST/functional/{tissue}_activities_pathways.csv',
        'results/ST/ST_{tissue}_deconvoluted.csv'
    output:
        view = 'results/ST/Misty/{tissue}/{sample}/pathwaysCT_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule run_views:
    input:
        view = 'results/ST/Misty/{tissue}/{sample}/{view_type}_view.rds'
    output: 
        directory('results/ST/Misty/{tissue}/{sample}/{view_type}_misty_model')
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
        lambda w: expand('results/ST/Misty/{{tissue}}/{sample}/{{view_type}}_misty_model', sample = config['samples'][w.tissue])
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

rule get_dif_interactions:
    input:
        'data/original/ST/metadata_visium_{tissue}.csv',
        lambda w: expand('results/ST/Misty/{{tissue}}/{sample}/{{view_type}}_misty_model', sample = config['samples'][w.tissue])
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
        view = 'results/ST/Misty/{tissue}/{sample}/{view_type}_view.rds'
    output:
        corr = temp('results/ST/Misty/{tissue}/{sample}/Corr_{view_type}.csv')
    params:
        corr = 'pearson'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/get_interactions_corr.R"

rule combine_interaction_corr:
    input:
        lambda w: expand('results/ST/Misty/{{tissue}}/{sample}/Corr_{{view_type}}.csv', sample = config['samples'][w.tissue])
    output:
        corr = 'results/ST/Misty/{tissue}/{view_type}_Corr.csv'
    shell:
        "awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} >> {output}"


def dif_interactions_inputs(wildcards):
    files = {'data': 'results/ST/{wildcards.tissue}_wImages.h5ad'.format(wildcards=wildcards),\
            'importances': 'results/Misty/{wildcards.tissue}/{wildcards.view_type}_importances.csv'.format(wildcards=wildcards),\
            'diffInteractions': 'results/Misty/{wildcards.tissue}/{wildcards.view_type}_diffInteractions.csv'.format(wildcards=wildcards),
            'correlations': 'results/ST/Misty/{wildcards.tissue}/{wildcards.view_type}_Corr.csv'.format(wildcards=wildcards)\
            }
    if (wildcards.view_type == 'celltype' or wildcards.view_type == 'pathwaysCT'):
        files['cellprops'] = 'results/ST/ST_{wildcards.tissue}_deconvoluted.csv'.format(wildcards=wildcards)

    if (wildcards.view_type == 'functional' or wildcards.view_type == 'pathwaysCT'):
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
