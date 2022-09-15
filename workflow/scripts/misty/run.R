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
  
  view_fp <- snakemake@input$view
  cores <- snakemake@threads[[1]]
  output_dir <- snakemake@output[[1]]
}else{
  #files for testing in Rstudio
  rs <- 42
  
  view_fp <- snakemake@input$view
  cores <- 6
  output_dir <- "mistyTest"
}


# load files --------------------------------------------------------------

cat("DEBUG: setting multisession plan\n")
#defining parallelisation
plan(multisession, workers = cores)

cat("DEBUG: reading misty view from", view_fp, "\n")
misty.views <- readRDS(view_fp)



# run misty ---------------------------------------------------------------

cat("DEBUG: started running misty with seed", rs,"\noutput dir is:", output_dir, "\n")
misty.views %>% run_misty(results.folder= output_dir, seed = rs, verbose = FALSE)
cat("INFO: finished building misty models; stored in", output_dir, "\n")



