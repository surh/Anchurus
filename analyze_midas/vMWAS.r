#!/usr/bin/env Rscript
library(MatrixEQTL)
library(ggplot2)
library(argparser)

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
  
  # Optional arguments
  p <- argparser::add_argument(p, "--nrows", help = "Number of rows to read at a time",
                               default = 5000, type = "numeric")
  p <- argparser::add_argument(p, "--outdir", help = "Output directory",
                               default = "out/", type = "character")
  
  # Read arguments
  args <- argparser::parse_args(p)
  
  return(args)
}



#############

# Arguments
# args <- list(snps = "snps.txt",
#              covariates = "covariates.txt",
#              phenotype = "phenotype.txt",
#              outdir = "out/",
#              nrows = 1000)
args <- process_arguments()

if(!dir.exists(args$outdir))
  dir.create(args$outdir)

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