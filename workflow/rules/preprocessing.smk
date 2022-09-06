# rule ST_seurat_to_h5ad:
#     input:
#         data = 'data/original/ST/ST_{tissue}_annotated.rds'
#     output:
#         h5ad = 'data/working/ST/ST_{tissue}_annotated.h5ad'
#     conda:
#         "../envs/preprocessing.yml"
#     script:
#         "../scripts/preprocessing/RDS_to_h5ad.R"

rule MO_seurat_to_h5ad:
    input:
        data = 'data/original/MO/MO_{tissue}_annotated.rds'
    output:
        h5ad = directory('data/working/MO/{tissue}')
    resources:
        mem_mb=30000
    conda:
        "../envs/preprocessing.yml"
    script:
        "../scripts/preprocessing/RDS_to_h5ad.R"

rule ST_seurat_to_h5ad:
    input:
        data = 'data/original/ST/ST_{tissue}_annotated.rds'
    output:
        h5ad = directory('data/working/ST/{tissue}')
    resources:
        mem_mb=30000
    conda:
        "../envs/preprocessing.yml"
    script:
        "../scripts/preprocessing/RDS_to_h5ad.R"

rule MO_h5ad_to_h5mu:
    input:
        rules.MO_seurat_to_h5ad.output
    output:
        muad = 'data/working/MO/{tissue}.h5mu'
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/preprocessing/MO_to_mudata.py"

rule ST_combine_to_h5ad:
    input:
        rules.ST_seurat_to_h5ad.output
    output:
        muad = 'data/working/ST/{tissue}.h5ad'
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/preprocessing/ST_to_adata.py"

rule combine_visium:
    input:
        'data/working/ST/{tissue}.h5ad',
        expand('data/original/ST/visium_data_{{tissue}}/{sample}/filtered_feature_bc_matrix.h5')
    output:
        'data/working/ST/{tissue}_wImages.h5ad'
    conda:
        "../envs/astromouse.yml"
    script:
        '../scripts/preprocessing/annData_w_images.py'
