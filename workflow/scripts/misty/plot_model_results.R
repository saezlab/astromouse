library(tidyr)
library(dplyr)
library(ggplot2)
library(mistyR)
library(factoextra)


extract_contrast_interactions <- function (misty.results.from, misty.results.to, 
                                           views = NULL, cutoff.from = 1, cutoff.to = 1, 
                                           trim = -Inf, trim.measure = c(
                                             "gain.R2", "multi.R2","intra.R2",
                                             "gain.RMSE", "multi.RMSE", "intra.RMSE"
                                           )){
  
  trim.measure.type <- match.arg(trim.measure)
  
  #check that both results collections are properly formatted
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results.from)), 
                          msg = "The first provided result list is malformed. Consider using collect_results().")
  assertthat::assert_that(("improvements.stats" %in% names(misty.results.from)), 
                          msg = "The provided result list is malformed. Consider using collect_results().")
  assertthat::assert_that(("importances.aggregated" %in% names(misty.results.to)), 
                          msg = "The second provided result list is malformed. Consider using collect_results().")
  assertthat::assert_that(("improvements.stats" %in% names(misty.results.to)), 
                          msg = "The provided result list is malformed. Consider using collect_results().")
  
  
  if (is.null(views)) {
    #check that both result collections contain the same views
    assertthat::assert_that(rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::pull(.data$view) %>% unique(),
                                                    misty.results.to$importances.aggregated %>% dplyr::pull(.data$view) %>% unique())),
                            msg = "The requested views do not exist in both result lists.")
    
    views <- misty.results.from$importances.aggregated %>% dplyr::pull(.data$view) %>% unique()
  } else {
    #check that all views are present in both result collections
    assertthat::assert_that(all(views %in% (misty.results.from$importances.aggregated %>% dplyr::pull(.data$view))) &
                              all(views %in% (misty.results.to$importances.aggregated %>% dplyr::pull(.data$view))),
                            msg = "The requested views do not exist in both result lists.")
    
  }
  
  #check that both collections contain exactly the same target-predictor interactions in each view
  assertthat::assert_that(all(views %>% purrr::map_lgl(function(current.view) {
    
    rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
                              dplyr::pull(.data$Predictor) %>% unique(),
                            misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
                              dplyr::pull(.data$Predictor) %>% unique())) &
      
      rlang::is_empty(setdiff(misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
                                dplyr::pull(.data$Target) %>% unique(),
                              misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view) %>%
                                dplyr::pull(.data$Target) %>% unique()))
    
  })), msg = "Incompatible predictors and targets.")
  
  
  inv <- sign((stringr::str_detect(trim.measure.type, "gain") |
                 stringr::str_detect(trim.measure.type, "RMSE", negate = TRUE)) - 0.5)
  
  
  targets <- misty.results.from$improvements.stats %>% dplyr::filter(.data$measure == trim.measure.type, inv * .data$mean >= inv * trim) %>% dplyr::pull(.data$target)
  
  interactions <- views %>% purrr::map_dfr(function(current.view) {
    
    from.view.wide <- misty.results.from$importances.aggregated %>% dplyr::filter(.data$view == current.view, .data$Target %in% targets) %>%
      tidyr::pivot_wider(names_from = "Target", values_from = "Importance", -c(.data$view, .data$nsamples))
    
    to.view.wide <- misty.results.to$importances.aggregated %>% dplyr::filter(.data$view == current.view, .data$Target %in% targets) %>%
      tidyr::pivot_wider(names_from = "Target", values_from = "Importance", -c(.data$view, .data$nsamples))
    
    
    
    
    mask <- ((from.view.wide %>% dplyr::select(-.data$Predictor)) < cutoff.from) & ((to.view.wide %>% dplyr::select(-.data$Predictor)) >= cutoff.to)
    
    
    if(sum(mask, na.rm = TRUE) > 0){ 
      # assertthat::assert_that(sum(mask, na.rm = TRUE) > 0, msg = paste0("All values are cut off while contrasting."))
      
      
      masked <- ((to.view.wide %>% tibble::column_to_rownames("Predictor")) * mask)
      
      
      masked %>% dplyr::slice(which(masked %>% rowSums(na.rm = TRUE) > 0)) %>%
        dplyr::select(which(masked %>% colSums(na.rm = TRUE) > 0)) %>% tibble::rownames_to_column("Predictor") %>%
        tidyr::pivot_longer(names_to = "Target", values_to = "Importance", -.data$Predictor) %>% dplyr::filter(Importance > 0) %>%
        dplyr::mutate(view = current.view)
    }else{
      mat = matrix(ncol = 0, nrow = 0)
      df=data.frame(mat)
    }
    
  })
  
  return(interactions)
}

reformat_samples <- function(misty.results, view){
  
  results <- lapply(misty.results, function(x){
    if('sample' %in% colnames(x)){
      if(view == 'celltype' | view == 'CTpathways'){
        x$sample <- x$sample %>% dirname() %>% basename()  
      }else{
        x$sample <- x$sample %>% basename()  
      }
      
    }
    return(x)
  })
  
  return(results)
}

# define input and outputs ------------------------------------------------

cat("DEBUG: defining inputs, outputs, and script parameters\n")

if(exists("snakemake")){
  
  tissue <- snakemake@wildcards$tissue
  view <- snakemake@wildcards$view_type
  
  plot_params <- snakemake@params[[1]]
  
  metadata_fp <- snakemake@input[[1]]
  result_folders <- unlist(snakemake@input[2:length(snakemake@input)])
  
}else{
  tissue <- 'brain'
  view <- 'celltype'
  
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

if(view == 'functional' | view == 'pathwaysCT' | view == 'CTpathways'){
  intra_name <- 'intra_act'
  cleaning <- TRUE
  
}

if(view == 'celltype'){
  intra_name <- 'intra'
  cleaning <- FALSE
  
}

if (view == 'celltype' | view == 'CTpathways'){
  
  result_folders <- lapply(result_folders, function(folder){
    list.dirs(folder, recursive = FALSE)
  }) %>% unlist()
  
}


results <- lapply(result_folders %>% collect_results(), function(x){
  if('view' %in% colnames(x)){
    x$view <- gsub('[[:digit:]]+', '', x$view)
    x$view <- gsub('\\.$', '', x$view)
  }
  return(x)
}) %>% reformat_samples(., view = view)

imp.signature <- extract_signature(results, type = "importance", trim = 2)

imp.signature.pca <- prcomp(imp.signature %>% select(-sample))

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

if(view != 'functional'){
  results %>% plot_improvement_stats()

  results %>% plot_improvement_stats("intra.R2")

  results %>% plot_view_contributions(trim = 1)

  results %>% plot_interaction_heatmap(intra_name, trim = plot_params$trim, cutoff = plot_params$cutoff, clean = cleaning)

  results %>% plot_interaction_heatmap('para', trim = plot_params$trim, cutoff = plot_params$cutoff, clean = cleaning)
}

fviz_pca_var(imp.signature.pca,
             col.var = "cos2", select.var = list(cos2 = 15), repel = TRUE,
             gradient.cols = c("#666666", "#377EB8", "#E41A1C"), col.circle = NA
) + theme_classic()

# if(exists("snakemake")) dev.off()


# by condition plots ------------------------------------------------------


grouped.results <- lapply(levels(metadata$condition), function(group){
  
  group.samples <- metadata %>% filter(condition == group)
  
  if(view == 'celltype' | view == 'CTpathways'){
    keep <- which(result_folders %>% dirname() %>% basename() %in% group.samples$sample)
  }else{
    keep <- which(result_folders %>% basename() %in% group.samples$sample)
  }
  group.results <- result_folders[keep] %>% collect_results()
  
  group.results <- group.results %>% reformat_samples(., view = view)
})

names(grouped.results) <- levels(metadata$condition)

contrast.interactions <- names(grouped.results) %>%  purrr::map_dfr(function(to.group){
  
  from.group <- names(grouped.results)[which(!grepl(to.group, names(grouped.results)))]
  
  # get interactions that are only in to.group, but not from.group
  # using default cutoff of 1 (on Importance), as 'being present'
  interactions <- extract_contrast_interactions(grouped.results[[from.group]], grouped.results[[to.group]]) %>% 
    tidyr::unite(col = 'Interaction', .data$view, .data$Predictor, .data$Target, sep = '_')
  
  # extract the importances per sample for these interactions and add metadata
  importances <- results$importances %>% tidyr::unite(col = 'Interaction', .data$view, .data$Predictor, .data$Target, remove = FALSE) %>% 
    dplyr::filter(.data$Interaction %in% (interactions %>% dplyr::pull(.data$Interaction))) %>% dplyr::select(-.data$Predictor, -.data$Target) %>% dplyr::left_join(metadata, by = 'sample')
  
  # t-test over conditions and do BH p.value adjustment
  stats <- importances %>% dplyr::group_by(.data$Interaction) %>% rstatix::t_test(data =., Importance ~ condition) %>% 
    dplyr::left_join(importances %>% dplyr::select(.data$Interaction, .data$view) %>% dplyr::distinct(), by = 'Interaction') %>% dplyr::group_by(.data$view) %>%
    rstatix::adjust_pvalue(method = "BH") %>% dplyr::select(.data$view, .data$Interaction, .data$statistic, .data$p, .data$p.adj) %>% dplyr::rename(t.value = .data$statistic, p.value = .data$p) %>% dplyr::ungroup()
  
  stats %>% plyr::adply(.margins = 1, function(x){
    x$Interaction <- base::gsub(paste('^', x$view, '_' , sep=''), '', x$Interaction) 
    return(x) 
  }) %>% tidyr::separate(col=.data$Interaction, into = c('Predictor', 'Target'), sep = '_') %>% dplyr::mutate(only.in = to.group)
  
})

views <- contrast.interactions %>% dplyr::pull(.data$view) %>% unique() %>% sort()

views %>% purrr::walk(function(current.view){
  
  long.data <- contrast.interactions %>% dplyr::filter(.data$view == current.view) %>%
    dplyr::mutate(sig = -log10(.data$p.adj) * sign(.data$t.value), is.sig.05 = .data$p.adj < 0.05, is.sig.1 = .data$p.adj < 0.1)
  
  inter.plot <- long.data %>% ggplot2::ggplot(ggplot2::aes(x = .data$Predictor, y = .data$Target)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$sig)) + ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) + ggplot2::coord_equal() + 
    ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "#8DA0CB", midpoint = 0) + ggplot2::labs(fill='-log10(p)') +
    geom_tile(data = long.data %>% dplyr::filter(.data$is.sig.05), aes(fill = .data$sig, color = is.sig.05), size = 1) + 
    scale_color_manual(guide = FALSE, values = c(`TRUE` = "black")) +
    ggplot2::ggtitle(paste('Condition specific interactions in', current.view))
  
  print(inter.plot)
  
})
if(exists("snakemake")) dev.off()



if(exists("snakemake")){
  output_filenames <- snakemake@output[2:length(snakemake@output)]
}else{
  output_filenames <- c('test1.pdf', 'test2.pdf')
}

mapply(function(results, output_fp){

  if(exists("snakemake")) pdf(output_fp)
  if(view != 'functional'){
    results %>% plot_improvement_stats()

    results %>% plot_improvement_stats("intra.R2")

    results %>% plot_view_contributions(trim = 1)

    results %>% plot_interaction_heatmap(intra_name, trim = plot_params$trim, cutoff = plot_params$cutoff, clean = cleaning)

    results %>% plot_interaction_heatmap('para', trim = plot_params$trim, cutoff = plot_params$cutoff, clean = cleaning)
  }
  if(exists("snakemake")) dev.off()

  return()

}, grouped.results, output_filenames)








