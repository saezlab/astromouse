import scanpy as sc
from os.path import join, normpath

input_dir = snakemake.input[0]

### Load RNA modalities
print('INFO: Loading RNA assays')
rna = sc.read_h5ad(join(input_dir,'RNA.h5ad'))
rna.layers['counts'] = rna.layers['logcounts'].copy()
del rna.layers['logcounts']
print(rna)

sct = sc.read_h5ad(join(input_dir,'SCT.h5ad'))

#change layer name
sct.layers['SCT'] = sct.layers['logcounts'].copy()
del sct.layers['logcounts']


sct.raw = rna

#change keys to lower case
for obsm in list(sct.obsm.keys()):
    sct.obsm[obsm.lower()] = sct.obsm[obsm].copy()
    del sct.obsm[obsm]

print(sct)

sct.write_h5ad(snakemake.output[0])
print('INFO: Combined adata written to', snakemake.output[0])
