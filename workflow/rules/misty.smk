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
        view = 'data/working/ST/Misty/{sample}/{tissue}_functional_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"

rule get_deconv_views:
    input:
        'data/working/ST/Misty/brain_coordinates.csv',
        'data/working/ST/ST_brain_deconvoluted.csv'
    output:
        view = 'data/working/ST/Misty/{sample}/brain_CT_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_views.R"
