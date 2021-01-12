#!/usr/bin/env Rscript

# (C) Copyright 2019-2021 Sur Herrera Paredes
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

#################### FUNCTIONS ##################
process_arguments <- function(){
  p <- arg_parser(paste("produce alignments from metagenomes"))

  # Positional arguments
  p <- add_argument(p, "info",
                    help = paste("snps_info.txt file"),
                    type = "character")

  p <- add_argument(p, "alleles",
                    help = paste("snps_alleles.txt file"),
                    type = "character")

  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Directory for output"),
                    default = "output/",
                    type = "character")
  # p <- add_argument(p, "--type",
  #                   help = paste("Either a single directory that contains SNV files for",
  #                                "one species ('single'), or a directory that contains",
  #                                "multiple species sub-directories ('multi')."),
  #                   default = "single",
  #                   type = "character")
  p <- add_argument(p, "--snvs",
                    help = paste("SNV effect type to use. Either 'all', 'ns' or 's'."),
                    default = 'all',
                    type = "character")

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)

  # Process arguments
  # if(args$type == 'multi'){
  #   cat("\tReading multiple input directories\n")
  #   args$midas_dir <- list.dirs(args$midas_dir, recursive = FALSE)
  # }
  # # Adding parameters that user cannot modify.
  # args$depth_thres <- 1
  # args$freq_thres <- 1
  # args$keep_last_codon <- TRUE
  # args$missing_as <- "gap"

  if(all(args$snvs != c('all', 'ns', 's'))){
    stop("ERROR: snvs must be in c('all', 'ns', 's)", call. = TRUE)
  }


  return(args)
}


#' Read data
#'
#' @param info_file
#' @param alleles_file
#' @param snvs
#'
#' @return
#' @export
#'
#' @examples
read_data <- function(info_file, alleles_file, snvs = "all"){
  info <- read_tsv(info_file,
                   col_types = cols(site_id = col_character(),
                                    ref_id = col_character(),
                                    ref_pos = col_number(),
                                    gene_id = col_character(),
                                    snp_effect = col_character(),
                                    locus_type = col_character()))


  # Select SNVs by snp_effect
  if(snvs == "ns"){
    info <- info %>%
      filter(snp_effect == "non-synonymous")
  }else if(snvs == "s"){
    info <- info %>%
      filter(snp_effect == "synonymous")
  }else if(snvs != "all"){
    stop("ERROR: snvs must be one of 'all', 'ns' or 's'",
         call. = TRUE)
  }

  # Get SNP IDs from genes
  snps <- info %>%
    filter(locus_type == "CDS") %>%
    select(site_id) %>%
    unlist()

  # Select only SNVs from genes
  info <- info %>%
    filter(site_id %in% snps) %>%
    select(site_id, gene_id)

  # Get alleles
  alleles <- read_tsv(alleles_file,
                      col_types = cols(site_id = col_character(),
                                       .default = col_character()))
  alleles <- alleles %>%
    filter(site_id %in% snps)

  # Create final data and return
  alleles <- alleles %>%
    left_join(info, by = "site_id")

  return(alleles)
}
##################################

args <- process_arguments()
# setwd("/home/sur/micropopgen/exp/2021/today2")
# args <- list(info = "test_in/snps_info.txt",
#              alleles =  "test_in/snps_alleles.txt",
#              snvs = 'all',
#              outdir = "output")

# Load rest of libraries
library(tidyverse)
library(seqinr)

# Read allele data
alleles <- read_data(info_file = args$info,
                     alleles_file = args$alleles,
                     snvs = args$snvs)

# Create output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

# Create alignnments
Res <- alleles %>%
  split(.$gene_id) %>%
  map(function(a, outdir){
    gene <- unique(a$gene_id)

    # Create aln object
    aln <- a %>%
      select(-site_id, -gene_id) %>%
      as.list %>%
      map(~ replace_na(.x, "-"))

    # Write
    filename <- file.path(outdir, paste0(gene, ".aln.fasta"))
    write.fasta(aln, names = names(aln), file.out = filename)
  }, outdir = args$outdir)
