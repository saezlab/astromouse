rule plot_pathways:
    input:
        'data/working/ST/{tissue}_wImages.h5ad'
    output:
        'plots/functional/{tissue}_pathways.pdf'
    params:
        normalisation = config['functional'].get("normalisation", 'log1p'), #or SCT
        top_genes = config['functional'].get("pathways_top_gene", 300)
    conda:
        "../envs/astromouse.yml"
    resources:
        mem_mb=40000
    script:
        '../scripts/functional/plots_pathways.py'

rule get_pathways:
    input:
        data = 'data/working/ST/{tissue}_wImages.h5ad'
    output:
        act = 'data/working/ST/functional/{tissue}_activities_{network}.csv' #network being either 'pathways' or 'TFs'
    params:
        normalisation = config['functional'].get("normalisation", 'log1p'), #or SCT
        top_genes = config['functional'].get("pathways_top_gene", 300),
        TF_conf = config['functional'].get("TF_confidence", 'ABC'),
        method = config['functional'].get("pathways_method", 'mlm')
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/functional/compute_activities.py"