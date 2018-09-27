library(MatrixEQTL)
library(ggplot2)
library(dplyr)

snps_file <- "snps.txt"
covariates_file <- "covariates.txt"
phenotype_file <- "phenotype.txt"


## Settings
useModel <- modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
output_file_name <- "results.txt"
pvOutputThreshold <- 1;
errorCovariance <- numeric();


snps <- SlicedData$new();
snps$fileDelimiter <- "\t";      # the TAB character
snps$fileOmitCharacters <- "NA"; # denote missing values;
snps$fileSkipRows <- 1;          # one row of column labels
snps$fileSkipColumns <- 1;       # one column of row labels
snps$fileSliceSize <- 2000;      # read file in slices of 2,000 rows
snps$LoadFile(snps_file);
snps

gene <- SlicedData$new();
gene$fileDelimiter <- "\t";      # the TAB character
gene$fileOmitCharacters <- "NA"; # denote missing values;
gene$fileSkipRows <- 1;          # one row of column labels
gene$fileSkipColumns <- 1;       # one column of row labels
gene$fileSliceSize <- 2000;      # read file in slices of 2,000 rows
gene$LoadFile(phenotype_file);
gene

cvrt <- SlicedData$new();
cvrt$fileDelimiter <- "\t";      # the TAB character
cvrt$fileOmitCharacters <- "NA"; # denote missing values;
cvrt$fileSkipRows <- 1;          # one row of column labels
cvrt$fileSkipColumns <- 1;       # one column of row labels
cvrt$LoadFile(covariates_file);
cvrt


mwa <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = 'qqplot',
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

plot(mwa)
