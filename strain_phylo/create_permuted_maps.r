# setwd("/cashew/users/sur/exp/fraserv/2021/today")
library(tidyverse)

args <- list(outdir = "output",
             n_maps = 100)

full_map <- read_tsv("healthy_gut_continent.tsv")
full_map

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}


set.seed(7235)
for( d in list.dirs("filtered_alleles/", full.names = T, recursive = F)){
  # d <- list.dirs("filtered_alleles/", full.names = T, recursive = F)[2]
  f <- file.path(d, "snps_alleles.txt")
  
  if(!file.exists(f))
    next()
  
  spec <- basename(d)
  cat(spec, "\n")
  
  alleles <- read_tsv(f, col_types = cols(site_id = col_character()))
  samples <- setdiff(colnames(alleles), "site_id")
  
  spec_map <- full_map %>%
    filter(sample %in% samples)
  
  # Prepare output dir
  spec_dir <- file.path(args$outdir, spec)
  dir.create(spec_dir)
  
  for(i in 1:args$n_maps){
    spec_map_i <- spec_map %>%
      mutate(group = sample(group))
    
    f_i <- file.path(spec_dir, paste0("map_", i, ".txt"))
    write_tsv(spec_map_i, f_i)
  }
}

