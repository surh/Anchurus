#!/usr/bin/env Rscript

# (C) Copyright 2020 Sur Herrera Paredes
#
# This file is part of Anchurus.
#
# Anchurus is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Anchurus is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public Licens

library(argparser)
library(HMVAR)

process_arguments <- function(){
  p <- arg_parser(paste("Reads midas data and makes quasi-phased",
                        "genotype calls"))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste("Directory with results of midas_merge.py",
                                 "snps"),
                    type = "character")
  p <- add_argument(p, "output",
                    help = paste("Path to output file"),
                    type = "character")
  
  # optional argumens
  p <- add_argument(p, "--min_depth",
                    help = paste("Min sequence depth."),
                    type = "integer",
                    default = 5)
  p <- add_argument(p, "--min_snv_prop",
                    help = paste("Minimum proportion of quasi-phaseable",
                                 "SNVs per sample"),
                    type = "integer",
                    default = 0.8)
  p <- add_argument(p, "--maf_thres",
                    help = paste("Minor allele frequency threshold for",
                                 "quasi-phasing SNVs"),
                    type = "integer",
                    default = 0.2)
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  return(args)
}

args <- process_arguments()

res <- qp_genotypes(midas_dir = args$midas_dir, min_depth = args$min_depth,
             min_snv_prop = args$min_snv_prop,
             maf_thres = args$maf_thres)
if(nrow(res) > 0){
  readr::write_tsv(res, args$output, na = ".")
}else{
  cat("No samples passed the thresholds")
}
