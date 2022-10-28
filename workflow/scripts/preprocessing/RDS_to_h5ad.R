library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

if(exists("snakemake")){
  input_fp = normalizePath(snakemake@input$data)
  output_fp <- snakemake@output[[1]]
}else{
  input_fp = normalizePath("data/original/ST/ST_brain_annotated.rds")
  output_fp <- "data/working/MO/MO_brain_annotated"
}

data <- readRDS(input_fp)

if(grepl(".h5ad$", output_fp, ignore.case = TRUE)){
  
  sce <- as.SingleCellExperiment(data)
  writeH5AD(data, file = output_fp)
  
}else if(grepl(".csv$", output_fp, ignore.case = TRUE)){
  
  if(exists("snakemake")){
    assay <- as.character(snakemake@params$assay)
    cellprop_cutoff = as.numeric(snakemake@params$cellprop_cutoff)
  }else{
    assay <- 'hvg2000'
  }
  
  chosen <- grep(paste('^',assay,'$', sep = ''), names(data@assays), value = TRUE)
  
  if(length(chosen) == 0) stop('The assay ', assay, ' does not exist in the Seurat object. Choose one of the following:\n', paste(names(data@assays), collapse = ' '))
  
  dec <- t(as.matrix(Seurat::GetAssayData(data, slot = 'data', assay = chosen)))
  # dec[dec < cellprop_cutoff] <- 0
  
  write.csv(dec, file = output_fp)
  
}else{
  
  if(dir.exists(output_fp)){
    unlink(output_fp, recursive=TRUE)
  }
  dir.create(output_fp)
  
  assays <- base::intersect(names(data@assays), c('RNA', 'ATAC', 'SCT', 'SCT_CC', 'chromvar'))
  
  lapply(assays, function (x){
    cat('INFO: converting assay:', as.character(x), '\n')
    temp <- data
    temp@active.assay <- x
    keep <- grepl(paste('^',paste(c(paste(paste0(x, c('_pca', '_harmony', '_umap'))), 'pca', 'harmony', 'umap'), collapse = '$|^'), '$', sep = ''), names(temp@reductions))
    temp@reductions[!keep] <- NULL
    
    cat('INFO: with reductions', names(temp@reductions), '\n')

    temp <- as.SingleCellExperiment(temp)

    writeH5AD(temp, file = paste(output_fp, .Platform$file.sep, gsub('_', '-', x), '.h5ad', sep = '' ))
    
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
