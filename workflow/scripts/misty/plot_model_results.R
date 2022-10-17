library(tidyr)
library(dplyr)
library(ggplot2)
library(mistyR)
library(factoextra)


# define input and outputs ------------------------------------------------

cat("DEBUG: defining inputs, outputs, and script parameters\n")

if(exists("snakemake")){
  
  tissue <- snakemake@wildcards$tissue
  view <- snakemake@wildcards$view_type
  
  plot_params <- snakemake@params[[1]]
  
  metadata_fp <- snakemake@input[[1]]
  result_folders <- unlist(snakemake@input[2:length(snakemake@input)])
  
}else{
  tissue <- 'heart'
  view <- 'functional'
  
  plot_params <- list(trim = 1, cutoff = 1)
  
  samples <- list.files(paste('data/original/ST/visium_data', tissue, sep = '_')) %>% sort()
  
  result_folders <- paste0(paste('results/ST/Misty', tissue, sep = .Platform$file.sep), '/', samples, '/', view, '_misty_model')
  
  #files for testing in Rstudio
  metadata_fp <- paste('data/original/ST/metadata_visium_', tissue,'.csv', sep = '')
  
}


# load data ---------------------------------------------------------------

result_folders <-  result_folders %>% sort()
samples <- result_folders %>% dirname() %>% basename()

metadata <- read.csv(metadata_fp)

if(tissue == 'brain'){
  metadata <- metadata %>% unite(., 'sample', sample_name, sub_array, sep = '_') %>% 
    rename('mouse' = 'mice_id.bio_origin') %>% select(-sample_no) %>% arrange(sample)
  
  metadata$condition <- as.factor(metadata$condition)
  levels(metadata$condition) <- c('Flight', 'Control')
  
}else if(tissue == 'heart'){
  metadata <- metadata %>% rename('sample' = 'Sample.Folder') %>% select(sample, condition, mouse) %>% arrange(sample)
  
  metadata$condition <- as.factor(metadata$condition)
  levels(metadata$condition) <- c('Flight', 'Control')
  
  metadata$sample <- gsub('-', '_', metadata$sample)
}

results <- lapply(result_folders %>% collect_results(), function(x){
  if('view' %in% colnames(x)){
    x$view <- gsub('[[:digit:]]+', '', x$view)
    x$view <- gsub('\\.$', '', x$view)
  }
  return(x)
})

imp.signature <- extract_signature(results, type = "importance", trim = 1)

imp.signature.pca <- prcomp(imp.signature %>% select(-sample))

if(view == 'functional'){
  intra_name <- 'intra_act'
  cleaning <- TRUE
  
}else if (view == 'celltype'){
  intra_name <- 'intra'
  cleaning <- FALSE
    
}

# plots -------------------------------------------------------------------

if(exists("snakemake")) pdf(snakemake@output[[1]])
ggplot(
  left_join(bind_cols(sample = metadata %>%  pull(sample), imp.signature.pca$x), metadata, by = 'sample'),
  aes(x = PC1, y = PC2)
) +
  geom_point(aes(color = as.factor(condition)), size = 1) +
  labs(color = "Condition") +
  theme_classic()



ggplot(
  left_join(bind_cols(sample = metadata %>%  pull(sample), imp.signature.pca$x), metadata, by = 'sample'),
  aes(x = PC1, y = PC2)
) +
  geom_point(aes(color = as.factor(mouse)), size = 1) +
  labs(color = "Mouse") +
  theme_classic()


results %>% plot_improvement_stats()

results %>% plot_improvement_stats("intra.R2")

results %>% plot_view_contributions(trim = 1)

results %>% plot_interaction_heatmap(intra_name, trim = plot_params$trim, cutoff = plot_params$cutoff, clean = cleaning)

results %>% plot_interaction_heatmap('para', trim = plot_params$trim, cutoff = cutoff, clean = cleaning)


fviz_pca_var(imp.signature.pca,
             col.var = "cos2", select.var = list(cos2 = 15), repel = TRUE,
             gradient.cols = c("#666666", "#377EB8", "#E41A1C"), col.circle = NA
) + theme_classic()

if(exists("snakemake")) dev.off()


# by condition plots ------------------------------------------------------


lapply(levels(metadata$condition), function(group){
  
  if(exists("snakemake")){
    if(group == 'Flight') pdf(snakemake@output[[2]])
    if(group == 'Control') pdf(snakemake@output[[3]])
  }
  
  group.samples <- metadata %>% filter(condition == group)
  keep <- which(result_folders %>% dirname() %>% basename() %in% group.samples$sample)
  
  group_folders <- result_folders[keep]
  
  results <- group_folders %>% collect_results()
  
  results %>% plot_improvement_stats()
  
  results %>% plot_improvement_stats("intra.R2")
  
  results %>% plot_view_contributions(trim = 1)
  
  results %>% plot_interaction_heatmap(intra_name, trim = plot_params$trim, cutoff = plot_params$cutoff, clean = cleaning)
  
  results %>% plot_interaction_heatmap('para', trim = plot_params$trim, cutoff = plot_params$cutoff, clean = cleaning)
  
  if(exists("snakemake")) dev.off()
  
  return()
  
})



