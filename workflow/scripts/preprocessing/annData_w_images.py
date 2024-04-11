import scanpy as sc
import muon as mu
import os
import pandas as pd
from glob import glob


# Define input and output files
if 'snakemake' in locals():
    tissue = snakemake.wildcards[0]
    filtered_adata = snakemake.input[0]
    sample_paths = snakemake.input[1:]
else:
    tissue = 'heart'
    filtered_adata = 'data/working/ST/{0}.h5ad'.format(tissue)
    sample_paths = [os.path.join(folder, 'filtered_feature_bc_matrix.h5') for folder in sorted(glob(os.path.join('data/original/ST/visium_data_{0}/'.format(tissue), '*')))]


# Read in filtered visium data (originally from Seurat)
adata = sc.read_h5ad(filtered_adata)
print(adata)
adata.obs.head()


# make the dataset to id mapping df (just to keep track which sample went into the seurat object in which order)
if tissue == 'brain':
    mapping_df = pd.DataFrame({'sample_name' : adata.obs['section.name'], 'dataset_id': adata.obs['slide.order']}).drop_duplicates().set_index('sample_name')
    mapping_df['dataset_id'] = mapping_df['dataset_id'].astype(int)
elif tissue == 'heart':
    mapping_df = pd.DataFrame(adata.obs.filter(['sample_id', 'mouse', 'run']),dtype=str).drop_duplicates()
    mapping_df['mouse'] = mapping_df['mouse'].str.split('_', expand=True)[1].astype(str).str[1:]
    mapping_df['dataset_id'] = pd.DataFrame(mapping_df.index.values, index = mapping_df.index.values)[0].str.split('_',expand = True)[1]
    mapping_df['sample_name'] = mapping_df[['run', 'sample_id', 'mouse']].agg('_'.join, axis = 1)
    mapping_df = mapping_df.filter(['dataset_id', 'sample_name'], axis = 1).sort_values('sample_name', axis = 0).set_index('sample_name')

# Define which metadata columns to keep
if tissue == 'brain':
    keep_columns = ['nCount_RNA', 'nFeature_RNA', 'nCount_SCT', 'nFeature_SCT', 'mouse', 'condition', 'seurat_clusters', 'SCT_snn_res.0.44',  'SCT_snn_res.0.5', 'annot']
    adata.obs = adata.obs.rename(columns={"bio_origin": "mouse", 'sample_condition':'condition'})
elif tissue == 'heart':
    keep_columns = ['nCount_RNA', 'nFeature_RNA', 'nCount_SCT', 'nFeature_SCT', 'sample_id', 'mouse', 'condition', 'seurat_clusters', 'SCT_snn_res.0.6', 'SCT_snn_res.0.5', 'annot', 'percent.mt', 'percent.hb', 'percent.pc']
    adata.obs['mouse'] = adata.obs['mouse'].str.split('_', expand=True)[1].astype(str)

# Filter metadata columns
adata.obs = adata.obs.filter(keep_columns, axis = 1)

# Had to rewrite read function because the data structure given by collaborators was changed from the default SpaceRanger output....
from pathlib import Path
from typing import Union, Optional
import json
import pandas as pd
from matplotlib.image import imread
from anndata import (
    AnnData
)

def custom_read_visium(
    path: Union[str, Path],
    genome: Optional[str] = None,
    *,
    count_file: str = "filtered_feature_bc_matrix.h5",
    library_id: str = None,
    load_images: Optional[bool] = True,
    source_image_path: Optional[Union[str, Path]] = None,
) -> AnnData:
    """\
    Read 10x-Genomics-formatted visum dataset.
    In addition to reading regular 10x output,
    this looks for the `spatial` folder and loads images,
    coordinates and scale factors.
    Based on the `Space Ranger output docs`_.
    See :func:`~scanpy.pl.spatial` for a compatible plotting function.
    .. _Space Ranger output docs: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/overview
    Parameters
    ----------
    path
        Path to directory for visium datafiles.
    genome
        Filter expression to genes within this genome.
    count_file
        Which file in the passed directory to use as the count file. Typically would be one of:
        'filtered_feature_bc_matrix.h5' or 'raw_feature_bc_matrix.h5'.
    library_id
        Identifier for the visium library. Can be modified when concatenating multiple adata objects.
    source_image_path
        Path to the high-resolution tissue image. Path will be included in
        `.uns["spatial"][library_id]["metadata"]["source_image_path"]`.
    Returns
    -------
    Annotated data matrix, where observations/cells are named by their
    barcode and variables/genes by gene name. Stores the following information:
    :attr:`~anndata.AnnData.X`
        The data matrix is stored
    :attr:`~anndata.AnnData.obs_names`
        Cell names
    :attr:`~anndata.AnnData.var_names`
        Gene names
    :attr:`~anndata.AnnData.var`\\ `['gene_ids']`
        Gene IDs
    :attr:`~anndata.AnnData.var`\\ `['feature_types']`
        Feature types
    :attr:`~anndata.AnnData.uns`\\ `['spatial']`
        Dict of spaceranger output files with 'library_id' as key
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['images']`
        Dict of images (`'hires'` and `'lowres'`)
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['scalefactors']`
        Scale factors for the spots
    :attr:`~anndata.AnnData.uns`\\ `['spatial'][library_id]['metadata']`
        Files metadata: 'chemistry_description', 'software_version', 'source_image_path'
    :attr:`~anndata.AnnData.obsm`\\ `['spatial']`
        Spatial spot coordinates, usable as `basis` by :func:`~scanpy.pl.embedding`.
    """
    path = Path(path)
    adata = sc.read_10x_h5(path / count_file, genome=genome)

    adata.uns["spatial"] = dict()

    from h5py import File

    with File(path / count_file, mode="r") as f:
        attrs = dict(f.attrs)
    if library_id is None:
        library_id = str(attrs.pop("library_ids")[0], "utf-8")

    adata.uns["spatial"][library_id] = dict()

    if load_images:

        folder_content = os.listdir(path)
        pos = [file for file in folder_content if 'tissue_positions_list' in file]
        scales = [file for file in folder_content if 'scalefactors' in file]
        image = [file for file in folder_content if 'tissue_hires' in file]

        files = dict(
            tissue_positions_file=path / pos[0],
            scalefactors_json_file=path / scales[0],
            hires_image=path / image[0],
            lowres_image=path / image[0],
        )

        # check if files exists, continue if images are missing
        for f in files.values():
            if not f.exists():
                if not any(x in str(f) for x in ["hires_image", "lowres_image"]):
                    # logg.warning(
                    #     f"You seem to be missing an image file.\n"
                    #     f"Could not find '{f}'."
                    # )
                # else:
                    raise OSError(f"Could not find '{f}'")

        adata.uns["spatial"][library_id]['images'] = dict()
        for res in ['hires', 'lowres']:
            try:
                adata.uns["spatial"][library_id]['images'][res] = imread(
                    str(files[f'{res}_image'])
                )
            except Exception:
                raise OSError(f"Could not find '{res}_image'")

        # read json scalefactors
        adata.uns["spatial"][library_id]['scalefactors'] = json.loads(
            files['scalefactors_json_file'].read_bytes()
        )

        adata.uns["spatial"][library_id]["metadata"] = {
            k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
            for k in ("chemistry_description", "software_version")
            if k in attrs
        }

        # read coordinates
        positions = pd.read_csv(files['tissue_positions_file'], header=None)
        positions.columns = [
            'barcode',
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm['spatial'] = adata.obs[
            ['pxl_row_in_fullres', 'pxl_col_in_fullres']
        ].to_numpy()
        adata.obs.drop(
            columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
            inplace=True,
        )

        # put image path in uns
        if source_image_path is not None:
            # get an absolute path
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(
                source_image_path
            )

    return adata

#Load each dataset individually
adatas = []
for sample in sample_paths:
    sample_name = os.path.basename(os.path.dirname(sample))
    print(sample_name)
    data = custom_read_visium(os.path.dirname(sample), library_id=sample_name)
    data.var_names_make_unique(join = '.')

    data.obs['labels'] = [ind + '_' + str(mapping_df.loc[sample_name, 'dataset_id']) for ind in data.obs.index.values]
    data.obs = data.obs.set_index('labels').rename_axis(index = None)

    adatas.append(data)

#Merge 
adata_spatial = sc.concat(
    adatas,
    join = 'outer',
    label="library_id",
    uns_merge="unique",
    keys=[
        k
        for d in [adatas[ii].uns["spatial"] for ii in range(len(adatas))]
        for k, v in d.items()
    ]#,
    # index_unique="-",
)

#Filter cells and features by those in the filtered object (from Seurat)
adata_spatial.raw = adata_spatial
adata_spatial = adata_spatial[adata.obs.index, list(set(adata.var.index.values).intersection(set(adata_spatial.var.index)))]
adata_spatial.obs = adata_spatial.obs.merge(adata.obs, left_index=True, right_index=True)

# Add transformed counts and embeddings
adata_spatial.layers['SCT'] = adata.X.copy()
for obsm in list(adata.obsm.keys()):
    adata_spatial.obsm[obsm] = adata.obsm[obsm].copy()
adata_spatial

# Write out
if 'snakemake' in locals():
    adata_spatial.write(snakemake.output[0])