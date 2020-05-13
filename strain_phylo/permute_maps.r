#!/usr/bin/env Rscript

# (C) Copyright 2019 Sur Herrera Paredes
#
# This file is part of HMVAR.
#
# HMVAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HMVAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HMVAR.  If not, see <http://www.gnu.org/licenses/>.

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
