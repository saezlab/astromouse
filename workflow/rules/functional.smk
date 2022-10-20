
def activities_inputs(wildcards):
    files = {'data': 'results/ST/{wildcards.tissue}_wImages.h5ad'.format(wildcards=wildcards)}
    if wildcards.network == 'GRNs':
        files['net'] = 'results/MO/celloracle/{wildcards.tissue}/GRNs/'.format(wildcards=wildcards)
    return files

rule plot_pathways:
    input:
        unpack(activities_inputs)
    params:
        lambda w: config['functional'][w.network]
    output:
        'plots/functional/{tissue}_{network}.pdf'
    conda:
        "../envs/astromouse.yml"
    resources:
        mem_mb=40000
    script:
        '../scripts/functional/plots_pathways.py'

rule get_PDactivities:
    input:
        unpack(activities_inputs)
    output:
        act = 'results/ST/functional/{tissue}_activities_{network}.csv' #network being either 'pathways' or 'TFs'
    params:
        lambda w: config['functional'][w.network]
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/functional/compute_activities.py"