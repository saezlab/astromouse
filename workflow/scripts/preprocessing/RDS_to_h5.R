# loading libraries -------------------------------------------------------

#for loading data
library(Seurat)
library(SeuratDisk)

#for converting seurat objects to h5ad
library(hdf5r)
library(dplyr)


# helper functions --------------------------------------------------------

seurat_write_h5 <- function(seurat_obj = NULL, file = NULL, assays = NULL, reductions_map= list('UMAP'='UMAP')){

  if(is.null(file)){
    stop('No such file or directory')
  }

  h5 <- hdf5r::H5File$new(filename = file, mode = 'w')


  tryCatch({

    seurat_to_h5(seurat_obj = seurat_obj, h5 = h5, assays = assays, reductions_map = reductions_map)

  },
  error = function(e) print(e),
  finally = {
    h5$close_all()
  })

}


seurat_to_h5 <- function(seurat_obj, h5, assays, reductions_map){

  #save the relevant assay data

  h5assays <- h5$create_group('assays')
  lapply(assays, function(assay){
    slot <- 'counts'
    if(assay == 'SCT' | assay == 'chromvar') slot <- 'data'

    mat <- GetAssayData(seurat_obj, slot = slot, assay = assay)

    # stopifnot("Class of assay is not dgCMatrix" = 'dgCMatrix' %in% class(mat))

    h5mat <- h5assays$create_group(assay)

    if('dgCMatrix' %in% class(mat)){

      h5mat[['values']] <- slot(object = mat, name = 'x')
      h5mat[['indices']] <- slot(object = mat, name = 'i')
      h5mat[['indptr']] <- slot(object = mat, name = 'p')
      h5mat[['dims']] <- rev(slot(object = mat, name = 'Dim'))


      if (!is.null(slot(object = mat, name = 'Dimnames')[[1]])) {
        h5mat[['var_names']] <- slot(object = mat, name = 'Dimnames')[[1]]
      } else {
        h5mat[['var_names']] <- rownames(mat)
      }


      h5mat[['obs_names']] <- slot(object = mat, name = 'Dimnames')[[2]]
      h5attr(h5mat, 'datatype') <- 'SparseMatrix'
      h5attr(h5mat, 'seurat_slot') <- slot
      print(h5mat)

    }else{
      h5mat_matrix <- h5mat$create_group('matrix')
      for(k in 1:(ncol(mat))){
        if(any(is.numeric(mat[,k]), is.integer(mat[,k]))){

          h5mat_matrix[[colnames(mat)[k]]] <- mat[,k]
          h5attr(h5mat_matrix[[colnames(mat)[k]]], 'origin_dtype') = 'number'

        }
      }
      h5mat[['var_names']] <- rownames(mat)
      h5mat[['obs_names']] <- colnames(mat)
      h5attr(h5mat, 'datatype') <- 'matrix'
      h5attr(h5mat, 'seurat_slot') <- slot
      print(h5mat)
    }
  })

  # --- save the cell annotations
  df_to_h5(df = slot(seurat_obj, name = 'meta.data'), h5 = h5, gr_name = 'obs')

  # --- save gene information
  if(!is.null(rownames(seurat_obj))) {
    df_to_h5(df = data.frame(row.names = rownames(seurat_obj)), h5 = h5, gr_name = 'var')
  }

  reductions_to_h5(seurat_obj, h5, reductions_map)

  print('INFO: Successfully converted Seurat object to h5')
  print('INFO: final h5 file is:\n')
  print(h5)

}


df_to_h5 <- function(df, h5, gr_name=NULL){
  h5df <- h5$create_group(gr_name)
  h5df[['index']] = rownames(df)
  if(ncol(df)>0){
    h5df[['colnames']] = colnames(df)
  }
  # factor to levels,character to levels,logical to levels
  for(k in names(df)){
    if(is.factor(df[[k]])){
      h5df[[k]]<- as.integer(df[[k]]) - 1L # for 0 begin
      h5df[[paste0(k,'_levels')]]<- levels(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'category'
    }
    if(is.character(df[[k]])){
      str_to_lvl <- factor(df[[k]])
      h5df[[k]]<- as.integer(str_to_lvl) - 1L
      h5df[[paste0(k,'_levels')]]<- levels(str_to_lvl)
      h5attr(h5df[[k]], 'origin_dtype') = 'string'
    }
    if(is.logical(df[[k]])){
      h5df[[k]] <- as.integer(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'bool'
    }
    if(any(is.numeric(df[[k]]),is.integer(df[[k]]))){
      h5df[[k]] <- df[[k]]
      h5attr(h5df[[k]], 'origin_dtype') = 'number'
    }
  }
}

reductions_to_h5 <- function(seurat_obj, h5, reductions_map){

  h5red <- h5$create_group('reductions')

  lapply(names(reductions_map), function(reduction){
    if(reduction %in% names(slot(seurat_obj, name = 'reductions'))){

      mat <- slot(slot(seurat_obj, name = 'reductions')[[reduction]], 'cell.embeddings')

      h5mat <- h5red$create_group(reductions_map[[reduction]])

      for(k in 1:(length(ncol(mat)) + 1)){
        if(any(is.numeric(mat[,k]), is.integer(mat[,k]))){

          h5mat[[colnames(mat)[k]]] <- mat[,k]
          h5attr(h5mat[[colnames(mat)[k]]], 'origin_dtype') = 'number'

        }
      }
    }

  })
}


# define input and output paths -------------------------------------------
if(exists("snakemake")){
  input_fp = normalizePath(snakemake@input$data)
  output_fp <- snakemake@output[[1]]

  modality <- snakemake@params[[1]]

#   assays <- snakemake@params[[1]]
#   reductions  <- snakemake@params[[2]]

}else{
  input_fp = normalizePath("data/original/ST/ST_brain_annotated.rds")
  input_fp = normalizePath('data/original/MO/MO_brain_annotated.RData')
  output_fp <- "test.h5"

  reductions  <- list( 'umap'= 'umap', 'pca' = 'pca', 'harmony' = 'harmony')
  assays <- c('SCT', 'RNA')

  modality <- 'MO'

}

if(modality=='ST'){

  reductions  <- list( 'umap'= 'umap', 'pca' = 'pca', 'harmony' = 'harmony')
  assays <- c('SCT', 'RNA')

}else if(modality=='MO'){

  reductions  <- list()
  assays <- c('RNA', 'ATAC', 'chromvar')

}


# load data ---------------------------------------------------------------

if(endsWith(input_fp, '.h5Seurat')) seurat_obj <- LoadH5Seurat(input_fp, assays = c('counts', 'data'))
if(endsWith(input_fp, '.rds')) seurat_obj <- readRDS(input_fp)
if(endsWith(input_fp, '.RData')) {
  load(input_fp)
  seurat_obj  <- seurat_multi
  rm(seurat_multi)
}


# dev ---------------------------------------------------------------------

assays <- intersect(assays, Assays(seurat_obj))

if(length(assays) == 0) {
  stop('None of the assays provided can be found in the Seurat object')
}


seurat_write_h5(seurat_obj, output_fp, assays = assays, reductions_map = reductions)



