rule brain_pathways:
    input:
        'data/working/ST/{tissue}_wImages.h5ad'
    output:
        'plots/functional/{tissue}_pathways.pdf'
    params:
        normalisation = 'log1p', #or SCT
        top_genes = 300
    conda:
        "../envs/astromouse.yml"
    resources:
        mem_mb=40000
    script:
        '../scripts/functional/plots_pathways.py'