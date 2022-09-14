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
        coords = 'data/working/ST/Misty/{tissue}_coordinates.csv',
        pathways ='data/working/ST/functional/{tissue}_activities_pathways.csv',
        TFs = 'data/working/ST/functional/{tissue}_activities_TFs.csv'
    output:
        view = 'data/working/ST/Misty/{sample}/{tissue}_functional_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_func_views.R"

rule get_deconv_views:
    input:
        coords = 'data/working/ST/Misty/brain_coordinates.csv',
        deconv ='data/working/ST/ST_brain_deconvoluted.csv'
    output:
        view = 'data/working/ST/Misty/{sample}/brain_CT_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_CT_views.R"
