#!/usr/bin/env Rscript
library(tidyverse)
library(HMVAR)
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste(""))
  
  # Positional arguments
  p <- add_argument(p, "alleles",
                    help = paste("site by sample alleles table. Must have
                                 column site_id"),
                    type = "character")
  p <- add_argument(p, "info",
                    help = paste("info file per site. Must have columns",
                                 "site_id, gene_id, ref_id, ref_pos,",
                                 "snp_effect, and substitution"),
                    type = "character")
  p <- add_argument(p, "map",
                    help = paste("sample metadata file. Must have column",
                    "sample"),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--outdir",
                     help = paste("Output directory"),
                     type = "character",
                     default = "output")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  
  return(args)
}

args <- process_arguments()
# args <- list(info = "alleles/snps_info.txt",
#              alleles = "alleles/snps_alleles.txt",
#              map = "example_map.txt",
#              outdir = "output")


cat("Reading data...\n")
alleles <- readr::read_tsv(args$alleles,
                           col_types = readr::cols(site_id = readr::col_character()))
info <- readr::read_tsv(args$info,
                        col_types = readr::cols(site_id = readr::col_character()))
meta <- read_tsv(args$map,
                 col_types = cols(sample = col_character()))

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}else{
  stop("ERROR: outdir already exists")
}

cat("Finding singletons...\n")
Sing <- find_singletons(alleles = alleles,
                        info = info %>%
                          dplyr::select(site_id, ref_id, ref_pos, gene_id,
                                        snp_effect, substitution),
                        meta = meta)

cat("Writing singletons...\n")
write_tsv(Sing, file.path(args$outdir, "singletons.txt"))

cat("Testing singletons...\n")
Res <- test_singleton_enrichment(Dat = Sing, info = info)

cat("Writing test results...\n")
write_tsv(Res, file.path(args$outdir, "singleton_tests.txt"))


