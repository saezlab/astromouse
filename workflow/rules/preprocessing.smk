from glob import glob

# extract individual MO assays from a seurat object
rule MO_seurat_to_h5ad:
    input:
        data = 'data/original/MO/MO_{tissue}_annotated.rds'
    output:
        h5ad = directory('results/MO/{tissue}')
    resources:
        mem_mb=30000
    conda:
        "../envs/preprocessing.yml"
    script:
        "../scripts/preprocessing/RDS_to_h5ad.R"

# creates a mudata object with three assays: rna, atac and chromvar
# rna has counts in .X, and SCT / SCT_CC normalised data in respective layers
# atac has only counts
# chromvar also
rule MO_h5ad_to_h5mu:
    input:
        rules.MO_seurat_to_h5ad.output
    output:
        muad = 'results/MO/{tissue}.h5mu'
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/preprocessing/MO_to_mudata.py"

# extract ST assays from seurat object
rule ST_seurat_to_h5ad:
    input:
        data = 'data/original/ST/ST_{tissue}_annotated.rds'
    output:
        h5ad = temp(directory('results/ST/convert/{tissue}_annotated'))
    resources:
        mem_mb=50000,
        disk_mb=20000
    conda:
        "../envs/preprocessing.yml"
    script:
        "../scripts/preprocessing/RDS_to_h5ad.R"

# combine assays (previously in Seurat objects) into Anndata
rule ST_combine_to_h5ad:
    input:
        rules.ST_seurat_to_h5ad.output
    output:
        ad = temp('results/ST/convert/{tissue}.h5ad')
    conda:
        "../envs/astromouse.yml"
    script:
        "../scripts/preprocessing/ST_to_adata.py"

# extract stereoscope deconvolution from Seurat to .csv
rule ST_extract_deconv:
    input:
        data = 'data/original/ST/ST_brain_deconvoluted.rds'
    output:
        csv = 'results/ST/ST_brain_deconvoluted.csv'
    params:
        assay = config['deconvolution'].get("assay", 'hvg2000'),
        cellprop_cutoff = config['deconvolution'].get('cellprop_cutoff')
    # resources:
    #     mem_mb=30000
    conda:
        "../envs/preprocessing.yml"
    script:
        "../scripts/preprocessing/RDS_to_h5ad.R"

# from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
# HTTP = HTTPRemoteProvider()

# # download spaceranger output for visium slides
# checkpoint visium_data:
#     input:
#         HTTP.remote("data.mendeley.com/api/datasets/fjxrcbh672/draft/files/ab629966-ac1a-4b9e-8371-a739a14a859b?a=69394d54-235c-436e-be60-520cd2899517", keep_local=False)
#     output:
#         directory('data/original/ST/visium_data_brain/')
#     shell:
#         "(test -d data/ || mkdir data) && "
#         "(test -d data/original/ || mkdir data/original) && "
#         "(test -d data/original/ST/ || mkdir data/original/ST) && "
#         "(test -d temp_vis/ || mkdir temp_vis) && "
#         "unzip {input} 'visium_data/*' -d temp_vis && "
#         "mv temp_vis/visium_data {output} && "
#         "rm -r temp_vis/"


# Combines the ST data into one AnnData object
# the .raw contains the original count data -> WITHOUT ANY FILTERING of spots or genes
# the .X matrix contains counts based on the provided Seurat object, with spots and genes filtered
# the SCT layer countains SCT transformed data using Seurat#s SCTtransformt method
def get_samples_from_tissue(wildcards):
    # download_output = checkpoints.visium_data.get(**wildcards).output[0]
    # folders = sorted(glob(os.path.join(download_output, '*')))
    folders = sorted(glob(os.path.join('data/original/ST/visium_data_brain/', '*')))
    sample_paths = [os.path.join(folder, 'filtered_feature_bc_matrix.h5') for folder in folders]
    return sample_paths

# add images and unfiltered counts to ST data in Anndata format
rule combine_visium:
    input:
        'results/ST/convert/{tissue}.h5ad',
        get_samples_from_tissue
    output:
        'results/ST/{tissue}_wImages.h5ad'
    conda:
        "../envs/astromouse.yml"
    resources:
        mem_mb=40000
    script:
        '../scripts/preprocessing/annData_w_images.py'
