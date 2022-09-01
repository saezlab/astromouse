library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

if(exists("snakemake")){
  input_fp = normalizePath(snakemake@input$data)
  output_fp <- snakemake@output[[1]]
}else{
  input_fp = normalizePath("data/original/ST/ST_brain_annotated.rds")
  output_fp <- "data/working/ST/ST_brain_annotated.h5ad"
}

data <- readRDS(input_fp)
data <- as.SingleCellExperiment(data)

writeH5AD(data, file = output_fp)

if('SCT' %in% names(data@assays)) data@assays$SCT@SCTModel.list <- list()

SeuratDisk::SaveH5Seurat(data, filename = gsub('.h5ad', '.h5Seurat', output_fp))

SeuratDisk::Convert(gsub('.h5ad', '.h5Seurat', output_fp), dest = "h5ad")

unlink(gsub('.h5ad', '.h5Seurat', output_fp))
