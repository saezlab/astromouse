import scanpy as sc
import muon as mu
from os.path import join, normpath

input_dir = snakemake.input[0]

### Load RNA modalities
print('INFO: Loading RNA assays')
rna = sc.read_h5ad(join(input_dir,'RNA.h5ad'))
del rna.layers['logcounts']
print(rna)

sct = sc.read_h5ad(join(input_dir,'SCT.h5ad'))
sct = sct[rna.obs.index,sct.var.index.sort_values()] #sort cells according to RNA assay, and order genes alphabetically
sct.layers['SCT'] = sct.X #storing SCT data in layer
del sct.layers['logcounts']
print(sct)

sctcc = sc.read_h5ad(join(input_dir,'SCT-CC.h5ad'))
sctcc = sctcc[sct.obs.index, sct.var.index] #sort according to sct assay
del sctcc.layers['logcounts']
print(sctcc)

#storing count data in X matrix
sct.X = rna.X[sct.obs.index, sct.var.index] #sort and subset features according to sct assay
sct.layers['SCT_CC'] = sctcc.layers['SCT_CC']

for obsm in list(sctcc.obsm.keys()):
    sct.obsm[obsm] = sctcc.obsm[obsm]

print('INFO: Loading ATAC assay')
atac = sc.read_h5ad(join(input_dir,'ATAC.h5ad'))
atac = atac[sct.obs.index, :]
del atac.layers['logcounts']

print('INFO: Loading chromvar assay')
chrv = sc.read_h5ad(join(input_dir,'chromvar.h5ad'))
chrv = chrv[sct.obs.index, :]

mdata = mu.MuData({"rna": sct, "atac": atac, "chromvar": chrv})

print(mdata)

mdata.write(snakemake.output[0])
print('INFO: Combined mudata written to', snakemake.output[0])