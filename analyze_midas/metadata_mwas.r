library(MatrixEQTL)
library(ggplot2)
library(dplyr)

meta_file <- "metadata_qin2012.txt"
snp_file <- "snps_freq.txt"
outdir <- "out/"
snp_thres <- 0.05
covariates <- c("age")
phenotype <- c("Diabetic")
phenotype_ref <- "N"
meta_id <- "ID"

# Read data
meta <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(meta)
snps <- read.table(snp_file, header = TRUE, sep = "\t", row.names = 1)
head(snps)

# Filter snps
p1 <- ggplot(data.frame(Freqs = rowMeans(snps)),
             aes(x = Freqs)) +
  geom_histogram(bins = 20) +
  scale_x_continuous(limits = c(0,1)) +
  AMOR::theme_blackbox
p1

snps <- snps[ p1$data$Freqs > snp_thres, ]

# Select metadata
row.names(meta) <- meta %>% pull(!!meta_id)
meta <- meta[ !is.na(meta %>% pull(!!phenotype)), ] 
meta <- meta[ colnames(snps), ]
selected_samples <- intersect(colnames(snps), row.names(meta))
meta <- meta[selected_samples,]

snps <- snps[,selected_samples]
snps <- snps[ rowMeans(snps) > 0.05, ]


pheno <- meta[,phenotype, drop = FALSE]
pheno[ , phenotype ] <- relevel(factor(pheno[ , phenotype ]),
                                ref = phenotype_ref)
pheno[ , phenotype ] <- as.numeric(pheno[ , phenotype ]) - 1
pheno <- pheno %>% t
covs <- meta %>% select(!!covariates) %>% t


# Write input files
dir.create(outdir)
outdir <- "test/"
pheno_file <- paste0(outdir,"/phenotypes.txt")
snps_file <- paste0(outdir, "/snps.txt")
covs_file <- paste0(outdir, "/covariates.txt")
write.table(pheno, pheno_file, sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
write.table(snps, snps_file, sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
write.table(covs, covs_file, sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)


## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Output file name
output_file_name = "test/out.txt"

# Only associations significant at this level will be saved
pvOutputThreshold = 1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 6000;      # read file in slices of 2,000 rows
snps$LoadFile(snps_file);
snps

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(pheno_file);
gene

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covs_file)>0){
  cvrt$LoadFile(covs_file);
}
cvrt

## Run the analysis

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

# unlink(output_file_name);
dat <- data.frame(freq = as.numeric(snps$FindRow("236")$row),
                  age = as.numeric(cvrt$FindRow("age")$row),
                  pheno = as.numeric(gene$FindRow("Diabetic")$row))
dat
summary(lm(pheno ~ freq + age, dat))
