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

library(argparser)
process_arguments <- function(){
  p <- arg_parser(paste("Code that runs RERconverge test on a binary phenotype"))
  
  # Positional arguments
  p <- add_argument(p, "indir",
                    help = paste("Directory containing one all trees from a given species",
                                 "or a list of subdirectories (one per species) with each",
                                 "one containing all tree species of the corresponding",
                                 "species."),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--type",
                     help = paste("Indicates if <indir> corresponds to a single species",
                                  "(single), or to a directory with subdirectories for",
                                  "multiples species (multi)."),
                     type = "character",
                     default = "single")
  p <- add_argument(p, "--map_file",
                    help = paste("Path to file with map. Should have 'ID' and 'group'",
                                 "columns."),
                    type = "character",
                    default = "map.txt")
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--focal_phenotype",
                    help = paste("Phenotypic value to be tested. Internally, tree tips",
                                 "with this value will be set to one, and to zero",
                                 "otherwise. If nothing passed, the script assumes that",
                                 "there are only two values in the Group column of the",
                                 "<--map_file>."),
                    default = NULL,
                    type = "character")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(!(args$type %in% c('single', 'multi'))){
    stop("ERROR: --type must be 'single' or 'multi'", call. = TRUE)
  }
  
  # Set other arguments
  args$scale <- TRUE
  args$weight <- TRUE
  
  return(args)
}

library(tidyverse)
library(RERconverge)