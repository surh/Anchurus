#!/usr/bin/env Rscript
library(MatrixEQTL)
library(ggplot2)

# Arguments
args <- list(snps = "snps.txt",
             covariates = "covariates.txt",
             phenotype = "phenotype.txt",
             outdir = "out/",
             nrows = 1000)

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