rule ST_seurat_to_h5ad:
    input:
        data = 'data/original/ST/ST_{tissue}_annotated.rds'
    output:
        h5ad = 'data/working/ST/ST_{tissue}_annotated.h5ad'
    conda:
        "../envs/preprocessing.yml"
    script:
        "../scripts/preprocessing/RDS_to_h5ad.R"

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
