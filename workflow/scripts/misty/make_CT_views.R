# MISTy
library(mistyR)
library(future)

# data manipulation
library(tidyverse)
library(purrr)
library(distances)

if(exists("snakemake")){
  tissue <- snakemake@wildcards[[1]]
  
  coord_fp = normalizePath(snakemake@input$coords)
  deconv_fp <- normalizePath(snakemake@input$deconv)
  sample <- basename(dirname(snakemake@output[[1]]))
}else{
  tissue <- 'brain'
  
  coord_fp <- normalizePath(paste('data/working/ST/Misty/', tissue, '_coordinates.csv', sep=""))
  deconv_fp <- normalizePath(paste('data/working/ST/ST_', tissue, '_deconvoluted.csv', sep=""))
  sample <- "Sample_158_A1"
}



# load data ---------------------------------------------------------------

datas_fp <- list(coord = coord_fp, CT = deconv_fp)

datas <- lapply(datas_fp, function(fp){
  data <- read_csv(fp)
  data <- data %>% column_to_rownames(var = colnames(data)[1])
})


# filter data for specified sample ----------------------------------------

datas$coord <- datas$coord %>% filter(library_id == sample)

datas[2:length(datas)] <- lapply(datas[2:length(datas)], function(data){
  data <- data[rownames(datas$coord),]
})




# determine distance between spots ----------------------------------------
# based on image coordinates instead of rows and columns

centroid <- round(colMeans(datas$coord %>% select(array_col, array_row)))
spots <- datas$coord %>% filter(array_col %in% (centroid[1]-18):(centroid[1]+18) &  array_row %in% (centroid[2]-18):(centroid[2]+18)) %>% 
  select(x,y)

#calculate the distance between selected spots
dM <- dist(spots, diag = F, upper = F) %>% unique() %>% sort()

#shift by one and subtract
d <- abs(dM[1:(length(dM)-1)] - dM[2:length(dM)])

#average of shortest distances
radius <- dM[1:which(d > 0. * dM[1])[1]] %>% mean() %>% round()


cat('The distance between spots is', as.character(radius), '\n')
cat('Using', as.character(2*radius), 'as l parameter in paraview creation\n')

# make views --------------------------------------------------------------

CT.view <- create_initial_view(datas$CT) %>% add_paraview(datas$coord %>% select(x,y), l = radius * 2)

if(exists("snakemake")){
  saveRDS(CT.view, snakemake@output[[1]])
}


