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
  p <- argparser::add_argument(p, "snps", help = paste0("SNPs file."),
                               type = "character")
  p <- argparser::add_argument(p, "covariates", help = paste0("File with covariates."),
                               type = "character")
  p <- argparser::add_argument(p, "phenotype", help = paste0("File with phenotype."),
                               type = "character")
  # p <- argparser::add_argument(p, "nsamples", help = paste0("Number of samples"),
  #                              type = "integer")
  
  # Optional arguments
  p <- argparser::add_argument(p, "--outfile", help = "File with results",
                               default = "association_results.txt", type = "character")
  p <- argparser::add_argument(p, "--maf", help = paste0("Minor allele frequency threshold"),
                               default = 0.05, type = "numeric")
  p <- argparser::add_argument(p, "--plot", help = "Plot p-values", default = FALSE,
                               flag = TRUE)
  
  # Read arguments
  args <- argparser::parse_args(p)
  
  return(args)
}



#############

# Arguments
args <- list(snps = "merged.snps/Veillonella_sp_62404/snps_freq.txt",
             covariates = "covariates.txt",
             phenotype = "phenotype.txt",
             outfile = "association_results.txt",
             maf = 0.05,
             plot = FALSE)
# args <- process_arguments()

# Read data
# col_types <- paste0(c('c', rep('n', args$nsamples)), collapse = "")
# phenotype <- read_tsv(args$phenotype)
# covariates <- read_tsv(args$covariates, col_types = col_types)
# snps <- read_tsv(args$snps, col_types = col_types, n_max = 10) # read it to get the colnames
# 
# col_types <- paste0(c('c', rep('n', 149)), collapse = "")
# snps <- read_tsv(args$snps, col_types = col_types)
# snps

auto_read_tsv <- function(file){
  # Get header
  header <- read_tsv(file, n_max = 10)
  col_types <- paste0(c('c', rep('n', ncol(header) - 1)), collapse = "")
  tab <- read_tsv(file, col_types = col_types)
  
  return(tab)
}

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

# Reformat covariate data
dat <- covariates %>%
  rbind(phenotype %>% 
          mutate(Covariate = Phenotype) %>%
          select(-Phenotype)) %>%
  gather(Sample, Value, -Covariate) %>%
  spread(Covariate, Value, fill = NA)



Res <- chunk_association(snps)
head(Res)
Res

var <- "1000051"
var <- "119378"
d <- snps %>% filter(SNP == var)
d <- cbind(dat, Frequency = as.numeric(d[1,-1]))
d
m1 <- lm(f1, d)
m1
summary(m1)

Res %>% filter(SNP == "1000051")

# Set params
useModel <- modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
output_file_name <- paste0(args$outdir, "/results.txt")
pvOutputThreshold <- 1;
errorCovariance <- numeric();

# Load data
snps <- SlicedData$new();
snps$fileDelimiter <- "\t";      # the TAB character
snps$fileOmitCharacters <- "NA"; # denote missing values;
snps$fileSkipRows <- 1;          # one row of column labels
snps$fileSkipColumns <- 1;       # one column of row labels
snps$fileSliceSize <- args$nrows;      # read file in slices of 2,000 rows
snps$LoadFile(args$snps);

pheno <- SlicedData$new();
pheno$fileDelimiter <- "\t";      # the TAB character
pheno$fileOmitCharacters <- "NA"; # denote missing values;
pheno$fileSkipRows <- 1;          # one row of column labels
pheno$fileSkipColumns <- 1;       # one column of row labels
pheno$fileSliceSize <- args$nrows      # read file in slices of 2,000 rows
pheno$LoadFile(args$phenotype);

cvrt <- SlicedData$new();
cvrt$fileDelimiter <- "\t";      # the TAB character
cvrt$fileOmitCharacters <- "NA"; # denote missing values;
cvrt$fileSkipRows <- 1;          # one row of column labels
cvrt$fileSkipColumns <- 1;       # one column of row labels
if(length(args$covariates) > 0){
  cvrt$LoadFile(args$covariates);
}

tests <- Matrix_eQTL_main(snps = snps,
                          gene = pheno,
                          cvrt = cvrt,
                          output_file_name = output_file_name,
                          pvOutputThreshold = pvOutputThreshold,
                          useModel = useModel, 
                          errorCovariance = errorCovariance, 
                          verbose = TRUE,
                          pvalue.hist = TRUE,
                          min.pv.by.genesnp = FALSE,
                          noFDRsaveMemory = FALSE)