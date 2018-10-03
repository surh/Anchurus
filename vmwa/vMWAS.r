#!/usr/bin/env Rscript
# library(MatrixEQTL)
library(ggplot2)
library(argparser)
library(readr)
library(dplyr)
library(tidyr)

#############
#' Process command line arguments
#'
#' @return List of arguments
#' 
#' @author Sur Herrera Paredes
#' 
#' @importFrom argparser arg_parser add_argument parse_args
#' @export
process_arguments <- function(){
  p <- argparser::arg_parser(paste0("Runs variant MWA"))
  
  # Positional arguments
  p <- argparser::add_argument(p, "snps", help = paste0("SNPs file. First column must ",
                                                        "be named 'site_id'."),
                               type = "character")
  p <- argparser::add_argument(p, "covariates", help = paste0("File with covariates. First ",
                                                              "column must be namer 'Covariate'."),
                               type = "character")
  p <- argparser::add_argument(p, "phenotype", help = paste0("File with phenotype. First ",
                                                             "column must be named 'Phenotype'."),
                               type = "character")

  # Optional arguments
  p <- argparser::add_argument(p, "--outfile", help = "File with results",
                               default = "association_results.txt", type = "character")
  p <- argparser::add_argument(p, "--maf", help = paste0("Minor allele frequency threshold"),
                               default = 0.05, type = "numeric")
  p <- argparser::add_argument(p, "--permutations", help = paste0("Whether to permute and ",
                                                                  "how may permutations."),
                               default = 0, type = "integer")
  p <- argparser::add_argument(p, "--plot", help = "Plot p-values", default = FALSE,
                               flag = TRUE)
  p <- argparser::add_argument(p, "--lib", help = paste0("Location of code"),
                               default = "~/micropopgen/src/Anchurus/vmwa/",
                               type = "character")
  p <- argparser::add_argument(p, "--seed", help = paste0("Seed for permutations."),
                               default = 5743, type = "integer")
  
  # Read arguments
  args <- argparser::parse_args(p)
  
  return(args)
}

#' Test all snps
#'
#' @param snps snp x sample table. First column must be 'SNP'
#' @param phenotype phenotype x sample table. First column must
#' be 'Phenotype'
#' @param covariate covariate x smaple table. First column must
#' be 'Covariate'
#' @param f1 Formula. Left-hand side must be called 'Frequency'
#'
#' @return
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate select
#' @importFrom tidyr gather spread
#' @examples
make_test <- function(snps, phenotype, covariate, f1){
  # Reformat covariate data
  dat <- covariates %>%
    rbind(phenotype %>% 
            dplyr::mutate(Covariate = Phenotype) %>%
            dplyr::select(-Phenotype)) %>%
    tidyr::gather(Sample, Value, -Covariate) %>%
    tidyr::spread(Covariate, Value, fill = NA)
  
  # Make test
  Res <- chunk_association(snps)
  
  return(Res)
}
#############

# Arguments
args <- list(snps = "merged.snps/Streptococcus_salivarius_58022/snps_freq.txt",
             covariates = "covariates.txt",
             phenotype = "phenotype.txt",
             outfile = "association_results.txt",
             maf = 0.05,
             permutations = 10,
             plot = FALSE,
             lib = "~/micropopgen/src/Anchurus/vmwa/",
             seed = 5743)
# args <- process_arguments()

# Source
source(paste0(args$lib, "/functions.r"))

# Read data
# col_types <- paste0(c('c', rep('n', args$nsamples)), collapse = "")
# phenotype <- read_tsv(args$phenotype)
# covariates <- read_tsv(args$covariates, col_types = col_types)
# snps <- read_tsv(args$snps, col_types = col_types, n_max = 10) # read it to get the colnames
# 
# col_types <- paste0(c('c', rep('n', 149)), collapse = "")
# snps <- read_tsv(args$snps, col_types = col_types)
# snps

# Read snps covariates and phenotype
snps <- auto_read_tsv(args$snps) %>% select(SNP = site_id, everything())
phenotype <- auto_read_tsv(args$phenotype)
covariates <- auto_read_tsv(args$covariates)

# Select samples
samples <- intersect(intersect(colnames(phenotype), colnames(snps)), colnames(covariates))
if(length(samples) == 0)
  stop("ERROR: No intersection between SNPs, covariates and phenotype")
snps <- snps %>% select(SNP,samples)
phenotype <- phenotype %>% select(Phenotype, samples)
covariates <- covariates %>% select(Covariate, samples)

# Filter by MAF
ii <- snps %>% select(-SNP) %>% rowMeans >= args$maf
snps <- snps[ ii, ]
snps
# ii <- snps %>% select(-SNP) %>% colSums > 0
# snps <- snps[ ,c(TRUE, ii) ]
# snps <- 
# 

# Make formula
f1 <- formula(paste(phenotype$Phenotype,
                    paste(c("Frequency", covariates$Covariate),
                          collapse = "+"),
                    sep = "~"))

Res <- make_test(snps, phenotype, covariate, f1)


