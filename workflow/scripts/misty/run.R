# MISTy
library(mistyR)
library(future)

# data manipulation
library(tidyverse)
library(purrr)
library(distances)

# plotting
library(ggplot2)

# define inputs and outputs,  and params ----------------------------------
cat("DEBUG: defining inputs, outputs, and script parameters\n")

if(exists("snakemake")){
  #files and parameters input by snakemake
  rs <- snakemake@params$seed
  
  bypass_intra <- snakemake@params$bypass_intra
  
  view_fp <- snakemake@input$view
  cores <- snakemake@threads[[1]]
  output_dir <- snakemake@output[[1]]
  
  view_type <- snakemake@wildcards$view_type
}else{
  #files for testing in Rstudio
  rs <- 42
  
  bypass_intra <- FALSE
  
  view_fp <- 'results/ST/Misty/brain/Sample_158_C1/celltype_view.rds'
  cores <- 6
  output_dir <- "mistyTest"
  
  view_type <- 'celltype'
}


# load files --------------------------------------------------------------

cat("DEBUG: setting multisession plan\n")
#defining parallelisation
plan(multisession, workers = cores)

cat("DEBUG: reading misty view from", view_fp, "\n")
misty.views <- readRDS(view_fp)

lapply(names(misty.views), function(view){
  if(view != 'misty.uniqueid'){
    colnames(misty.views[[view]]$data) <<- gsub("[[:punct:]]", "", colnames(misty.views[[view]]$data))
  }
  return()
})


# run misty ---------------------------------------------------------------

cat("DEBUG: started running misty with seed", rs,"\nbypass.intra set to:", bypass_intra, "\noutput dir is:", output_dir, "\n")

if(view_type == 'celltype'){
  
  lapply(colnames(misty.views$intraview$data), function(target){
    
    temp.views <- misty.views %>% filter_views(NA, view = 'mask', .data[[target]]) %>% remove_views('mask')
    
    views <- names(misty.views)[!grepl("misty.uniqueid|mask", names(misty.views))]
    
    
    names(temp.views[views]) %>% lapply(function(current.view){
      
      keep.columns <- temp.views[[current.view]]$data %>% summarise(across(where(is.numeric), var)) %>% t() %>% as.data.frame() %>% 
        filter(V1 > 0) %>% row.names()
      
      cat(length(keep.columns), '\n')
      
      if(current.view == 'intraview'){
        keep.columns <- unique(c(keep.columns, target))
      }
      
      temp.views[[current.view]]$data <- temp.views[[current.view]]$data[keep.columns]
      
      print(ncol(temp.views[[current.view]]$data))
      
    })
    
    
    
    temp.views %>% run_misty(results.folder = paste(output_dir, target, sep = '/'), seed = rs, bypass.intra = bypass_intra, verbose = FALSE, target.subset = target)
  })
  
  
}else{
  
  misty.views %>% run_misty(results.folder = output_dir, seed = rs, bypass.intra = bypass_intra, verbose = FALSE)
}



cat("INFO: finished building misty models; stored in", output_dir, "\n")



