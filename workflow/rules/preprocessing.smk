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
        mem_mb=30000,
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
        ad = 'results/ST/convert/{tissue}.h5ad'
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
        assay = config['deconvolution'].get("assay", 'hvg2000')
    # resources:
    #     mem_mb=30000
    conda:
        "../envs/preprocessing.yml"
    script:
        "../scripts/preprocessing/RDS_to_h5ad.R"


# Combines the ST data into one AnnData object
# the .raw contains the original count data -> WITHOUT ANY FILTERING of spots or genes
# the .X matrix contains counts based on the provided Seurat object, with spots and genes filtered
# the SCT layer countains SCT transformed data using Seurat#s SCTtransformt method
def samples_from_tissue(wildcards):
    folders = sorted(glob(os.path.join('data/original/ST/visium_data_{0}/'.format(wildcards.tissue), '*')))
    sample_paths = [os.path.join(folder, 'filtered_feature_bc_matrix.h5') for folder in folders]
    return sample_paths

# add images and unfiltered counts to ST data in Anndata format
rule combine_visium:
    input:
        'results/ST/convert/{tissue}.h5ad',
        samples_from_tissue
    output:
        'results/ST/{tissue}_wImages.h5ad'
    conda:
        "../envs/astromouse.yml"
    resources:
        mem_mb=40000
    script:
        '../scripts/preprocessing/annData_w_images.py'
