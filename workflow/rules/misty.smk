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
        'results/ST/functional/{tissue}_activities_pathways.csv',
        'results/ST/functional/{tissue}_activities_TFs.csv'
    params:
        skip = 'intra'
    output:
        view = 'results/ST/Misty/{tissue}/{sample}/functional_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_deconv_views:
    input:
        'results/ST/Misty/brain_coordinates.csv',
        'results/ST/ST_brain_deconvoluted.csv'
    output:
        view = 'results/ST/Misty/brain/{sample}/celltype_view.rds'
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
        "../scripts/misty/test.R"

rule plot_misty_results:
    input:
        'data/original/ST/metadata_visium_brain.csv',
        lambda w: expand('results/ST/Misty/{{tissue}}/{sample}/{{view_type}}_misty_model', sample = config['samples'][w.tissue])
    output: 
        'plots/Misty/{tissue}/{view_type}_misty.pdf'
    params:
        seed = config['misty'].get("random_seed", 42)
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/test.R"

