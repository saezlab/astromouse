library(tidyr)
library(dplyr)
library(ggplot2)
library(mistyR)
library(factoextra)


# define input and outputs ------------------------------------------------

cat("DEBUG: defining inputs, outputs, and script parameters\n")

if(exists("snakemake")){
  
  tissue <- print(snakemake@wildcards$tissue)
  
  metadata_fp <- snakemake@input[[1]]
  result_folders <- unlist(snakemake@input)
  
}else{
  tissue <- 'brain'
  
  #files for testing in Rstudio
  metadata_fp <- 'data/original/ST/metadata_visium_brain.csv'
  result_folders <- c("data/working/ST/Misty/brain/Sample_159_B1/celltype_misty_model", "data/working/ST/Misty/brain/Sample_159_A1/celltype_misty_model")
  
}


# load data ---------------------------------------------------------------

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
}

results <- lapply(result_folders %>% sort() %>% collect_results(), function(x){
  if('view' %in% colnames(x)){
    x$view <- gsub('[[:digit:]]+', '', x$view)
    x$view <- gsub('\\.$', '', x$view)
  }
  return(x)
})

imp.signature <- extract_signature(results, type = "importance", trim = 1)

imp.signature.pca <- prcomp(imp.signature %>% select(-sample))


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

results %>% plot_interaction_heatmap('intra', trim = 1)

results %>% plot_interaction_heatmap('para', trim = 1)


fviz_pca_var(imp.signature.pca,
             col.var = "cos2", select.var = list(cos2 = 15), repel = TRUE,
             gradient.cols = c("#666666", "#377EB8", "#E41A1C"), col.circle = NA
) + theme_classic()

if(exists("snakemake")) dev.off()







