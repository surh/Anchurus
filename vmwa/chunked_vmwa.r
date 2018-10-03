#!/usr/bin/env Rscript

# (C) Copyright 2018 Sur Herrera Paredes
# 
# This file is part of Anchurus
# 
# Anchurus is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Anchurus is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Anchurus.  If not, see <http://www.gnu.org/licenses/>.

library(ggplot2)
# library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(argparser)
##########################################################
#' Internal function
#'
#' @param d 
#' @param f1 
#'
#' @return
#'
#' @examples
fit_model <- function(d, f1) {
  #snp <- unique(d$SNP)
  #cat(snp, "\n")
  m1 <- tryCatch(lm(f1, data = d, singular.ok = FALSE),
                 error = function(e){list(coefficients=NA)} )
  if( is.na(m1$coefficients[1]) || is.na(coef(m1)['Frequency']) ){
    r <- rep(NA, 4)
    names(r) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  }else{
    r <- summary(m1)$coefficients["Frequency", ]
  }
  
  return(r)
}

#' Chunk association
#' 
#' Function to be used with readr chunk readers. It requires global variables
#' dat and f1. NEED TO MAKE MORE MODULAR.
#'
#' @param snps A table of SNP frequencies per table.
#' @param pos Current chunk position
#'
#' @return A data.frame with results of linear regression
#' 
#' @author Sur Herrera Paredes
#'
#' @examples
chunk_association <- function(snps, pos){
  # Reformat snps
  snps <- snps %>% gather(Sample, Frequency, -SNP, na.rm = TRUE)
  # cat(dim(snps), "\n")
  
  # Merge snps with covariates
  dat <- snps %>% inner_join(dat, by = "Sample")
  # cat(dim(dat), "\n")
  
  # Apply linear regression
  res <- plyr::ddply(dat, "SNP", fit_model, f1 = f1)
  return(res)
}

#' p-value qqplot
#' 
#' quantile-quantile plot of expected vs observed
#' p-values.
#'
#' @param pvector Vector of p-values
#' @param main Title
#' @param ... 
#'
#' @return
#' @export
#' 
#' @references 
#' \url{http://GettingGeneticsDone.blogspot.com/}
#' \url{http://gettinggeneticsdone.blogspot.com/p/copyright.html}
#'
#' @examples
ggd.qqplot = function(pvector, main=NULL, ...) {
  # http://GettingGeneticsDone.blogspot.com/
  # See http://gettinggeneticsdone.blogspot.com/p/copyright.html
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

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
  p <- argparser::add_argument(p, "nsamples", help = paste0("Number of samples"),
                               type = "integer")
  
  # Optional arguments
  p <- argparser::add_argument(p, "--chunk_size", help = "Number of rows to read at a time",
                               default = 2000, type = "numeric")
  p <- argparser::add_argument(p, "--outfile", help = "File with results",
                               default = "association_results.txt", type = "character")
  p <- argparser::add_argument(p, "--plot", help = "Plot p-values", default = FALSE,
                               flag = TRUE)
  
  # Read arguments
  args <- argparser::parse_args(p)
  
  return(args)
}
########################################################

#setwd("/godot/users/sur/exp/fraserv/2018/today/")
#args <- list(snps = "test.txt",
#              covariates = "covariates.txt",
#              phenotype = "phenotype.txt",
#              chunk_size = 2000,
#              nsamples = 368,
#              outfile = "association_results.txt",
#              plot = TRUE)
args <- process_arguments()

# Read data
col_types <- paste0(c('c', rep('n', args$nsamples)), collapse = "")
phenotype <- read_tsv(args$phenotype)
covariates <- read_tsv(args$covariates, col_types = col_types)
snps <- read_tsv(args$snps, col_types = col_types, n_max = 10) # read it to get the colnames
#snps


# Check name consistency
if(any(colnames(snps)[-1] != colnames(covariates)[-1]))
  stop("ERROR: Covariates do not match samples")
if(any(colnames(snps)[-1] != colnames(phenotype)[-1]))
  stop("ERROR: Phenotypes do not match samples")

### Preprocessing
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

# Process file by chunk
#snps <- read_tsv(args$snps, col_types = col_types)
#snps

#snps <- snps %>% gather(Sample, Frequency, -SNP, na.rm = TRUE)
#snps
#unique(snps$SNP)
# cat(dim(snps), "\n")
  
# Merge snps with covariates
#dat <- snps %>% inner_join(dat, by = "Sample")
#dat
#unique(dat$SNP)
# cat(dim(dat), "\n")
  
# Apply linear regression
#res <- plyr::ddply(dat, "SNP", fit_model, f1 = f1)

#d <- dat %>% filter(SNP == "G437.327699")

#fit_model(d, f1)


Res <- read_tsv_chunked(file = args$snps,
                        callback = DataFrameCallback$new(chunk_association),
                        chunk_size = args$chunk_size,
                        col_types = col_types) %>% select(SNP, beta = Estimate,
                                                          SE = "Std. Error",
                                                          t.value = "t value",
                                                          p.value = "Pr(>|t|)")

# Write output
write_tsv(Res, args$outfile)

if(args$plot){
  library(ggplot2)
  library(svglite)
  
  p1 <- ggplot(Res, aes(x = p.value)) +
    geom_histogram(bins = 20) +
    theme(panel.background = element_rect(color = "black", fill = NA, size= 3),
          panel.grid = element_blank())
  ggsave("p.value_histogram.svg", p1, width = 6, height = 4)
  
  svglite::svglite("p.value_qqplot.svg", width = 6, height = 6)
  ggd.qqplot(Res$p.value)
  dev.off()
}
