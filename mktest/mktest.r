#!/usr/bin/env Rscript
library(tidyverse)
library(HMVAR)
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Perform McDonald-Kreitman test on allele tables"))
  
  # Positional arguments
  p <- add_argument(p, "alleles",
                    help = paste("site by sample allele table"),
                    type = "character")
  p <- add_argument(p, "info",
                    help = paste("SNV info table"),
                    type = "character")
  p <- add_argument(p, "map",
                    help = paste("mapping file."),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--output",
                     help = paste("File to write results"),
                     type = "character",
                     default = "mkres.txt")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

args <- process_arguments()

# Read data
cat("reading data...\n")
info <- readr::read_tsv(file = args$info,
                        col_types = readr::cols(site_id = readr::col_character(),
                                                locus_type = readr::col_character(),
                                                gene_id = readr::col_character(),
                                                ref_id = readr::col_character()))
alleles <- readr::read_tsv(args$alleles,
                           col_types = readr::cols(site_id = readr::col_character()))
meta <- readr::read_tsv(args$map,
                        col_types = readr::cols(sample = readr::col_character(),
                                                group = readr::col_character()))

# Run mktest
cat("running mktest...\n")
Res <- mktest(alleles = alleles, info = info, map = meta)
cat("writing mktest...\n")
write_tsv(Res, args$output)
