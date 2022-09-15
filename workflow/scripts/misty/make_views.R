# MISTy
library(mistyR)
library(future)

# data manipulation
library(tidyverse)
library(purrr)
library(distances)

if(exists("snakemake")){
  tissue <- snakemake@wildcards[[1]]
  
  datas_fp <- lapply(snakemake@input, normalizePath)
  
  cat(snakemake@input)
  
  cat('\n\n\n\n')
  
  cat(datas_fp, '\n\n\n')
  
  if(length(datas_fp) == 3){
    # names(datas_fp) <- c('coord', 'pathways', 'TFs')
    
    view <- 'activities'
    
  }else if(length(datas_fp) == 2){
    # names(datas_fp) <- c('coord', 'CT')
    
    view <- 'cell_types'
    
  }else{
    stop('The input has to contain either 2 or 3 files, corresponding (in order) to coordinates, and then either intra+para data (in 1 file), or intra and paraview (2 files). ', length(datas_fp), ' files were given:\n', paste(datas_fp, collapse = ' '))
  }
  
  sample <- basename(dirname(snakemake@output[[1]]))
  
}else{
  tissue <- 'heart'
  
  coord_fp <- normalizePath(paste('data/working/ST/Misty/', tissue, '_coordinates.csv', sep=""))
  pathways_fp <- normalizePath(paste('data/working/ST/functional/', tissue, '_activities_pathways.csv', sep=""))
  tf_fp <- normalizePath(paste('data/working/ST/functional/', tissue, '_activities_TFs.csv', sep=""))
  sample <- "V19T26_014_A1_G8"
  view <- 'activities'
  
  datas_fp <- list(coord = coord_fp, pathways = pathways_fp, TFs = tf_fp)
  
}



# load data ---------------------------------------------------------------



datas <- lapply(datas_fp, function(fp){
  data <- read_csv(fp)
  data <- data %>% column_to_rownames(var = colnames(data)[1])
})


# filter data for specified sample ----------------------------------------

datas$coord <- datas$coord %>% filter(library_id == sample)

datas[2:length(datas)] <- lapply(datas[2:length(datas)], function(data){
  data <- data[rownames(datas$coord),]
})

if(nrow(datas[[1]]) < 1) stop('There are no spots for sample ', sample, ' in the provided coordinates file:\n', coord_fp)


# determine distance between spots ----------------------------------------
# based on image coordinates instead of rows and columns

centroid <- round(colMeans(datas$coord %>% select(array_col, array_row)))
spots <- datas$coord %>% filter(array_col %in% (centroid[1]-18):(centroid[1]+18) &  array_row %in% (centroid[2]-18):(centroid[2]+18)) %>% 
                                  select(x,y)

if(nrow(datas[[1]]) < 30) spots <- datas$coord

#calculate the distance between selected spots
dM <- dist(spots, diag = F, upper = F) %>% unique() %>% sort()

#shift by one and subtract
d <- abs(dM[1:(length(dM)-1)] - dM[2:length(dM)])

#average of shortest distances
radius <- dM[1:which(d > 0.1 * dM[1])[1]] %>% mean() %>% round()


cat('The distance between spots is', as.character(radius), '\n')
cat('Using', as.character(2*radius), 'as l parameter in paraview creation\n')

# make views --------------------------------------------------------------

if(view == 'activities'){
  intra.view <- create_initial_view(datas[[2]])
  
  para.view <- create_initial_view(datas[[3]]) %>% add_paraview(datas$coord %>% select(x,y), l = radius * 2)
  para.view <- within(para.view, rm(misty.uniqueid, intraview))
  
  misty.views <- intra.view %>% add_views(new.views = para.view)
  
}else if (view == 'cell_types'){
  
  misty.views <- create_initial_view(datas[[2]]) %>% add_paraview(datas$coord %>% select(x,y), l = radius * 2)
  
}


if(exists("snakemake")){
  saveRDS(misty.views, snakemake@output[[1]])
}


