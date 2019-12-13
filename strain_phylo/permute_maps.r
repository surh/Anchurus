#!/usr/bin/env Rscript
# setwd("/cashew/users/sur/exp/fraserv/2019/today2")

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Takes a mapping file and permutes the group ccolumn"))
  
  # Positional arguments
  p <- add_argument(p, "map_file",
                    help = paste("Path of the mapping file. A tab-delimited",
                                 "file with columns ID and Group."),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                     help = paste("Path to save output. Will fail if already exists."),
                     type = "character",
                     default = "output/")
  p <- add_argument(p, "--nperms",
                    help = "Number of permutations",
                    default = 100,
                    type = "integer")
  p <- add_argument(p, "--seed",
                    help = paste("Seed for randomization"),
                    type = "integer",
                    default = NULL)
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}



library(tidyverse)
args <- process_arguments()
# args <- list(map_file = "../2019-10-30.spec_maps/gut.maps/Bacteroides_vulgatus_57955.map.txt",
#              nperms = 3,
#              outdir = "output",
#              seed = 237)

map <- read_tsv(args$map_file)  

if(dir.exists(args$outdir)){
  stop("ERROR: output directory already exists.", call. = TRUE)
}else{
  dir.create(args$outdir)
}

if(!is.na(args$seed)){
  set.seed(args$seed)
}
Res <- 1:args$nperms %>%
  purrr::map(~write_tsv(x = tibble(ID = map$ID, Group = sample(map$Group, replace = FALSE)),
                        path = file.path(args$outdir, paste0("map_", ., ".txt"))))
  

