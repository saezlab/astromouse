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
  
  print(snakemake@input)
  
  cat('\n\n\n\n')
  
  print(datas_fp)
  
  if(length(datas_fp) == 3){
    
  }else if(length(datas_fp) == 2){
    
  }else{
    stop('The input has to contain either 2 or 3 files, corresponding (in order) to coordinates, and then either intra+para data (in 1 file), or intra and paraview (2 files). ', length(datas_fp), ' files were given:\n', paste(datas_fp, collapse = ' '))
  }
  
  sample <- snakemake@wildcards$sample
  
  if('cellprop_cutoff' %in% names(snakemake@params)){
    cellprop_cutoff <- snakemake@params$cellprop_cutoff
  }
  
  output_fp <- snakemake@output[[1]]
  
}else{
  tissue <- 'brain'
  view_type <- 'CTpathways'
  
  coord_fp <- normalizePath(paste('results/ST/Misty/', tissue, '_coordinates.csv', sep=""))
  celltype_fp <- normalizePath(paste('results/ST/ST_', tissue, '_deconvoluted.csv', sep=""))
  pathways_fp <- normalizePath(paste('results/ST/functional/', tissue, '_activities_pathways.csv', sep=""))
  tf_fp <- normalizePath(paste('results/ST/functional/', tissue, '_activities_TFs.csv', sep=""))
  sample <- "Sample_159_B1"
  
  output_fp <- paste('results/ST/Misty/', tissue,'/', view_type, '/views/', sample,'_view.rds', sep = '')
  
  cellprop_cutoff <- 0.05
  
  datas_fp <- list(coord = coord_fp, cellprop = celltype_fp)
  datas_fp <- list(coord = coord_fp, cellprop = celltype_fp, pathways = pathways_fp)
  
}


# determine view type -----------------------------------------------------

view <- output_fp %>% dirname() %>% dirname() %>% basename()



# load data ---------------------------------------------------------------



datas <- lapply(datas_fp, function(fp){
  data <- read_csv(fp)
  data <- data %>% column_to_rownames(var = colnames(data)[1])
})


# filter data for specified sample ----------------------------------------

datas[[1]] <- datas[[1]] %>% filter(library_id == sample)

datas[2:length(datas)] <- lapply(datas[2:length(datas)], function(data){
  data <- data[rownames(datas[[1]]),]
})

if(nrow(datas[[1]]) < 1) stop('There are no spots for sample ', sample, ' in the provided coordinates file:\n', datas_fp[[1]])


# determine distance between spots ----------------------------------------
# based on image coordinates instead of rows and columns

centroid <- round(colMeans(datas[[1]] %>% select(array_col, array_row)))

#%>% filter(array_col %in% (centroid[1]-18):(centroid[1]+18) &  array_row %in% (centroid[2]-18):(centroid[2]+18)) %>% 
spots <- datas[[1]] %>% select(x,y)

#calculate the distance between selected spots
dM <- dist(spots, method = "euclidean", diag = TRUE, upper = TRUE, p = 2) %>% unique() %>% sort()

#shift by one and subtract
d <- abs(dM[1:(length(dM)-1)] - dM[2:length(dM)])

#average of shortest distances
radius <- dM[1:which(d > 0.1 * dM[1])[1]] %>% mean() %>% round()


cat('The distance between spots is', as.character(radius), '\n')
cat('Using', as.character(2*radius), 'as l parameter in paraview creation\n')

# make views --------------------------------------------------------------

# TODO: redo this for all three views for this project
if(view == 'functional' | view == 'pathwaysCT' | view == 'CTpathways'){
  intra.view <- create_initial_view(datas[[2]])
  
  para.view <- create_initial_view(datas[[3]]) %>% add_paraview(datas[[1]] %>% select(x,y), l = radius * 2)
  para.view <- para.view %>% rename_view(., old.name = 'intraview', new.name = 'intra_act')
  para.view <- within(para.view, rm(misty.uniqueid))
  
  misty.views <- intra.view %>% add_views(new.views = para.view)
  
}

if (view == 'celltype'| view == 'CTpathways'){
  
  if (view == 'celltype'){
    misty.views <- create_initial_view(datas[[2]]) %>% add_paraview(datas[[1]] %>% select(x,y), l = radius * 2)
  }
  
  mask <- misty.views$intraview$data %>% transmute_all(function(x){x >= cellprop_cutoff})
  
  misty.views <- misty.views %>% add_views(., create_view('mask', mask))
  
}

para.name <- names(misty.views)[grepl('para', names(misty.views))]

misty.views <- misty.views %>% rename_view(., old.name = para.name, new.name = 'paraview', new.abbrev = 'para')

para <- misty.views$paraview$data
rownames(para) <- rownames(datas[[1]])

if(exists("snakemake")){
  saveRDS(misty.views, snakemake@output[[1]])
  write.csv(para, snakemake@output[[2]])
}

