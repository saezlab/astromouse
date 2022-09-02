library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

if(exists("snakemake")){
  input_fp = normalizePath(snakemake@input$data)
  output_fp <- snakemake@output[[1]]
}else{
  input_fp = normalizePath("data/original/MO/MO_brain_annotated.rds")
  output_fp <- "data/working/MO/MO_brain_annotated.h5ad"
}

data <- readRDS(input_fp)

if(grepl(".h5ad$", output_fp, ignore.case = TRUE)){
  
  sce <- as.SingleCellExperiment(data)
  writeH5AD(data, file = output_fp)
  
}else{
  
  if(dir.exists(output_fp)){
    unlink(output_fp, recursive=TRUE)
  }
  dir.create(output_fp)
  
  assays <- c('RNA', 'ATAC', 'SCT', 'SCT_CC', 'chromvar')
  
  lapply(assays, function (x){
    print('INFO: converting assay:', x)
    temp <- data
    temp@active.assay <- x
    keep <- grepl(paste('^',paste(paste0(x, c('_pca', '_harmony', '_umap')), collapse = '$|^'), '$', sep = ''), names(temp@reductions))
    temp@reductions[!keep] <- NULL
    
    print('INFO: with reductions', temp@reductions)

    sce <- as.SingleCellExperiment(temp)
    
    writeH5AD(data, file = paste(output_fp, .Platform$file.sep, gsub('_', '-', x), '.h5ad', sep = '' ))
    
  })
  
}




# 
# if('SCT' %in% names(data@assays)) data@assays$SCT@SCTModel.list <- list()
# 
# SeuratDisk::SaveH5Seurat(data, filename = gsub('.h5ad', '.h5Seurat', output_fp))
# 
# SeuratDisk::Convert(gsub('.h5ad', '.h5Seurat', output_fp), dest = "h5ad")
# 
# unlink(gsub('.h5ad', '.h5Seurat', output_fp))
