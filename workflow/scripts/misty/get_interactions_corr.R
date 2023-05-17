library(mistyR)
library(tibble)
library(tidyr)
library(dplyr)
library(purrr)
library(stringr)


# files -------------------------------------------------------------------

if(exists("snakemake")){
  #files and parameters input by snakemake
  
  view_type <- snakemake@wildcards$view_type
  sample <- snakemake@wildcards$sample
  tissue <- snakemake@wildcards$tissue
  
  correlation.type <- snakemake@params$corr
  
  interactions_fp <- snakemake@input$interactions
  view_fp <- snakemake@input$view

  output_fp <- snakemake@output[[1]]
  
  
}else{
  #files for testing in Rstudio
  
  view_type <- 'celltype'
  sample <- 'Sample_158_A1'
  tissue <- 'brain'
  
  correlation.type <- 'pearson'
  
  interactions_fp <- str_glue('results/Misty/{tissue}/{view_type}_diffInteractions.csv', tissue = tissue, view_type=view_type, .sep = "")
  view_fp <- str_glue('results/ST/Misty/{tissue}/{sample}/{view_type}_view.rds', tissue = tissue, sample=sample, view_type=view_type, .sep = "")
  
  output_fp <- 'test.csv'
}


# load data ---------------------------------------------------------------

interactions <- read.csv(interactions_fp) %>% select(view, Target, Predictor, model)
misty.views <- readRDS(view_fp)



# do correlation spot by spot ---------------------------------------------

targets <- intersect(interactions$Target, colnames(misty.views$intraview$data))


output <- targets %>%   map_dfr(function(target){
  
  current.inter <- interactions %>% filter(.data$Target == target)
  
  if(view_type == 'celltype' | view_type == 'CTpathways'){
    temp.views <- misty.views %>% filter_views(NA, view = 'mask', .data[[target]]) %>% remove_views('mask') 
  }else{
    temp.views <- misty.views
  }

  if(view_type == 'CTpathways'){
    print(colnames(temp.views$intra_act$data))
    colnames(temp.views$intra_act$data)  <- gsub("-", "", colnames(temp.views$intra_act$data))
    colnames(temp.views$paraview$data)  <- gsub("-", "", colnames(temp.views$paraview$data))
    print(colnames(temp.views$intra_act$data))
  }
  
  inter_cor <- current.inter %>% pull(.data$view) %>% unique() %>% map_dfr(function(current.view){
    
    #select columns of target-predictors
    if(current.view == 'intra'){
      df <- temp.views$intraview$data
    }else if (current.view == 'intra_act'){
      df <- cbind(temp.views$intraview$data %>% select(all_of(target)),
                  temp.views$intra_act$data %>% select(any_of(current.inter %>%
                                                                filter(.data$view == current.view) %>% pull(.data$Predictor) %>% unique())))
    }else if (current.view == 'para'){
      df <- cbind(temp.views$intraview$data %>% select(all_of(target)),
                  temp.views$paraview$data %>% select(any_of(current.inter %>%
                                                                filter(.data$view == current.view) %>% pull(.data$Predictor) %>% unique())))
    }
  
    correlations <- df %>%
      cor(method = correlation.type) %>%
      as.data.frame %>%
      rownames_to_column(var = 'var1') %>%
      gather(var2, value, -var1) %>% rename(corr = value)
    
    inner_join(current.inter %>% filter(.data$view == current.view), correlations, by = c('Predictor' = 'var1', 'Target' = 'var2'))
    
  })
  
  
  
  
}) %>% mutate(sample = sample)

if(exists("snakemake")){
  write.csv(output, file = output_fp, sep = ",")
}
