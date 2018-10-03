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
  p <- argparser::add_argument(p, "--plot", help = "Plot p-values",
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
  Res <- association(snps, dat = dat, f1 = f1)
  
  return(Res)
}


association <- function(snps, dat = dat, f1 = f1){
  # Reformat snps
  snps <- snps %>% tidyr::gather(Sample, Frequency, -SNP, na.rm = TRUE)
  # cat(dim(snps), "\n")
  
  # Merge snps with covariates
  dat <- snps %>% dplyr::inner_join(dat, by = "Sample")
  # cat(dim(dat), "\n")
  
  # Apply linear regression
  res <- plyr::ddply(dat, "SNP", fit_model, f1 = f1)
  return(res)
}
#############

# Arguments
# args <- list(snps = "merged.snps/Streptococcus_salivarius_58022/snps_freq.txt",
#              covariates = "covariates.txt",
#              phenotype = "phenotype.txt",
#              outfile = "association_results.txt",
#              maf = 0.05,
#              permutations = 10,
#              plot = FALSE,
#              lib = "~/micropopgen/src/Anchurus/vmwa/",
#              seed = 5743)
args <- process_arguments()

# Source
source(paste0(args$lib, "/functions.r"))

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

# Make formula
f1 <- formula(paste(phenotype$Phenotype,
                    paste(c("Frequency", covariates$Covariate),
                          collapse = "+"),
                    sep = "~"))

Res <- make_test(snps = snps, phenotype = phenotype, covariate = covariate, f1 = f1)

if(args$permutations > 0){
  Res$N <- 1
  Res$P <- 1
  set.seed(seed = args$seed)
  
  for(i in 1:args$permutations){
    if( i %% 100 == 1 ) cat("Running permutation ", i, "...\n", sep = '')
    
    # Permute phenotype
    pheno.p <- phenotype
    colnames(pheno.p) <- c(colnames(pheno.p)[1],sample(colnames(pheno.p)[-1], replace = FALSE))
    pheno.p <- pheno.p %>% select(Phenotype, samples)
    
    # Generate null stat
    res.p <- make_test(snps = snps, phenotype = pheno.p, covariate = covariates, f1 = f1)
    
    # Count
    cmp <- abs(res.p$beta) >= abs(Res$beta)
    na_count <- !is.na(cmp)
    Res$N <- Res$N + na_count
    cmp[!na_count] <- 0
    Res$P <- Res$P + cmp
  }
  
  Res$P <- Res$P / Res$N
}

# Write results
write_tsv(Res, args$outfile)



if(args$plot){
  library(ggplot2)
  # library(svglite)
  
  p1 <- ggplot(Res, aes(x = p.value)) +
    geom_histogram(bins = 20) +
    theme(panel.background = element_rect(color = "black", fill = NA, size= 3),
          panel.grid = element_blank())
  ggsave("p.value_histogram.svg", p1, width = 6, height = 4)
  
  png("p.value_qqplot.png", width = 1000, height = 1000, res = 300)
  ggd.qqplot(Res$p.value)
  dev.off()
  
  if(args$permutations > 0){
    p1 <- ggplot(Res, aes(x = P)) +
      geom_histogram(bins = 20) +
      theme(panel.background = element_rect(color = "black", fill = NA, size= 3),
            panel.grid = element_blank())
    ggsave("P_histogram.svg", p1, width = 6, height = 4)
    
    png("P_qqplot.png", width = 1000, height = 1000, res = 300)
    ggd.qqplot(Res$P)
    dev.off()
  }
}
