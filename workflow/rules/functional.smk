rule brain_pathways:
    input:
        'data/working/ST/brain_wImages.h5ad'
    output:
        'plots/functional/brain_pathways.pdf'
    params:
        normalisation = 'log1p' #or SCT
        top_genes = 300
    conda:
        "../envs/astromouse.yml"
    # resources:
    #     mem_mb=40000
    script:
        '../scripts/preprocessing/plots_pathways_brain.py'