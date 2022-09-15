rule get_coords:
    input:
        data = 'data/working/ST/{tissue}_wImages.h5ad'
    output:
        coords = 'data/working/ST/Misty/{tissue}_coordinates.csv'
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/misty/extract_coordinates.py"

rule get_func_views:
    input:
        'data/working/ST/Misty/{tissue}_coordinates.csv',
        'data/working/ST/functional/{tissue}_activities_pathways.csv',
        'data/working/ST/functional/{tissue}_activities_TFs.csv'
    output:
        view = 'data/working/ST/Misty/{tissue}/{sample}/functional_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_deconv_views:
    input:
        'data/working/ST/Misty/brain_coordinates.csv',
        'data/working/ST/ST_brain_deconvoluted.csv'
    output:
        view = 'data/working/ST/Misty/brain/{sample}/celltype_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule run_func_views:
    input:
        view = 'data/working/ST/Misty/{tissue}/{sample}/{view_type}_view.rds'
    output: 
        directory('data/working/ST/Misty/{tissue}/{sample}/{view_type}_misty_model')
    params:
        seed = config['misty'].get("random_seed", 42)
    conda:
        "../envs/misty.yml"
    threads: 6
    resources:
        mem_mb=25000,
        disk_mb=1000,
        time='12:00:00'
    script:
        "../scripts/misty/make_views.R"
