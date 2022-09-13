rule get_coords:
    input:
        data = 'data/working/ST/{tissue}_wImages.h5ad'
    output:
        coords = 'data/working/ST/Misty/{tissue}_coordinates.csv'
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/misty/extract_coordinates.py"

rule get_views:
    input:
        coords = 'data/working/ST/Misty/{tissue}_coordinates.csv',
        pathways ='data/working/ST/functional/{tissue}_activities_pathways.csv',
        TFs = 'data/working/ST/functional/{tissue}_activities_TFs.csv'
    output:
        view = 'data/working/ST/Misty/{sample}/{tissue}_misty_view.rds'
    conda:
        "../envs/misty.yml"
    script:
        "../scripts/misty/make_func_views.R"

