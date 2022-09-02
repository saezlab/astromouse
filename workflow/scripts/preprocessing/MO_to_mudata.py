import scanpy as sc
import muon as mu
from os.path import join, normpath

input_dir = snakemake.input[0]

### Load RNA modalities
print('INFO: Loading RNA assays')
rna = sc.read_h5ad(join(input_dir,'RNA.h5ad'))
rna.layers['counts'] = rna.layers['logcounts'].copy()
del rna.layers['logcounts']
print(rna)

sct = sc.read_h5ad(join(input_dir,'SCT.h5ad'))
sct.layers['SCT'] = sct.layers['logcounts'].copy()
del sct.layers['logcounts']
print(sct)

sctcc = sc.read_h5ad(join(input_dir,'SCT-CC.h5ad'))
sctcc.layers['SCT_CC'] = sctcc.layers['logcounts'].copy()
del sctcc.layers['logcounts']
print(sctcc)

sct.raw = rna
sct.layers['SCT_CC'] = sctcc.layers['SCT_CC']

for obsm in list(sctcc.obsm.keys()):
    sct.obsm[obsm] = sctcc.obsm[obsm]

print('INFO: Loading ATAC assay')
atac = sc.read_h5ad(join(input_dir,'ATAC.h5ad'))
atac.layers['counts'] = atac.layers['logcounts'].copy()
del atac.layers['logcounts']


print('INFO: Loading chromvar assay')
chrv = sc.read_h5ad(join(input_dir,'chromvar.h5ad'))


mdata = mu.MuData({"rna": sct, "atac": atac, "chromvar": chrv})

print(mdata)

mdata.write(snakemake.output[0])
print('INFO: Combined mudata written to', snakemake.output[0])