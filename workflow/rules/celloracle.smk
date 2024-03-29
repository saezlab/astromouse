rule peak_corr:
    input:
        data="results/MO/{tissue}.h5mu"
    singularity:
        "workflow/envs/celloracle.sif"
    output:
        path_all_peaks="results/MO/celloracle/{tissue}/all_peaks.csv",
        path_connections="results/MO/celloracle/{tissue}/cicero_connections.csv",
        path_plot="results/MO/celloracle/{tissue}/peak_thr.pdf"
    params:
        organism=config['organism'],
        min_count=config['celloracle']['min_count'],
        max_count=config['celloracle']['max_count']
    resources:
        mem_mb=48000
    shell:
        """
        Rscript workflow/scripts/celloracle/peak_corr.R {input.data} {params.organism} {params.min_count} {params.max_count} {output.path_plot} {output.path_all_peaks} {output.path_connections}
        """

rule tss_annotation:
    input:
        all_peaks="results/MO/celloracle/{tissue}/all_peaks.csv",
        connections="results/MO/celloracle/{tissue}/cicero_connections.csv"
    singularity:
        "workflow/envs/celloracle.sif"
    output:
        "results/MO/celloracle/{tissue}/processed_peak_file.csv"
    params:
        organism=config['organism'],
        thr_coaccess=config['celloracle']['thr_coaccess']
    resources:
        mem_mb=48000
    shell:
         "python workflow/scripts/celloracle/tss_annotation.py -a {input.all_peaks} -c {input.connections} -o {params.organism} -t {params.thr_coaccess} -p {output}"

rule tf_motif_scan:
    input:
        "results/MO/celloracle/{tissue}/processed_peak_file.csv"
    singularity:
        "workflow/envs/celloracle.sif"
    output:
        "results/MO/celloracle/{tissue}/motifs.celloracle.tfinfo"
    resources:
        mem_mb=48000
    params:
        organism=config['organism'],
        fpr=config['celloracle']['fpr']
    shell:
        "python workflow/scripts/celloracle/tf_motif_scan.py -p {input} -o {params.organism} -f {params.fpr} -t {output}"

rule build_base_grn:
    input:
        "results/MO/celloracle/{tissue}/motifs.celloracle.tfinfo"
    singularity:
        "workflow/envs/celloracle.sif"
    params:
        thr_motif_score=config['celloracle']['thr_motif_score']
    output:
        "results/MO/celloracle/{tissue}/base_GRN_dataframe.csv"
    resources:
        mem_mb=48000
    shell:
        "python workflow/scripts/celloracle/build_base_grn.py -i {input} -t {params.thr_motif_score} -g {output}"

rule build_grn:
    input:
        mdata="results/MO/{tissue}.h5mu",
        base_grn="results/MO/celloracle/{tissue}/base_GRN_dataframe.csv"
    singularity:
        "workflow/envs/celloracle.sif"
    output:
        "results/MO/celloracle/{tissue}/grn.celloracle.links"
    # resources:
    #     mem_mb=48000
    shell:
        "python workflow/scripts/celloracle/build_grn.py -m {input.mdata} -b {input.base_grn} -t {wildcards.tissue} -l {output}"

rule filter_grn:
    input:
        "results/MO/celloracle/{tissue}/grn.celloracle.links"
    singularity:
        "workflow/envs/celloracle.sif"
    params:
        thr_edge_pval=config['celloracle']['thr_edge_pval'],
        thr_top_edges=config['celloracle']['thr_top_edges'],
    output:
        path_out= directory('results/MO/celloracle/{tissue}/GRNs/')
    # resources:
    #     mem_mb=48000
    shell:
        """
        python workflow/scripts/celloracle/filter_grn.py -l {input} -p {params.thr_edge_pval} -t {params.thr_top_edges} -o {output.path_out}
        """

