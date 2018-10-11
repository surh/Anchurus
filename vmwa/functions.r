library(readr)

#' Automatic reading of tsv
#' 
#' Reads a SNP frequency, covariate or phenotype files.
#' Assumes first column is an ID column and the rest of
#' columns contain numeric values. Automatically figures
#' out number of samples.
#'
#' @param file File path or connection
#'
#' @return A tibble
#' @export
#' 
#' @importFrom readr read_tsv
#' @examples
auto_read_tsv <- function(file){
  # Get header
  header <- readr::read_tsv(file, n_max = 10)
  col_types <- paste0(c('c', rep('n', ncol(header) - 1)), collapse = "")
  tab <- readr::read_tsv(file, col_types = col_types)
  
  return(tab)
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
#' @importFrom plyr ddply
#' @importFrom dplyr inner_join
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @examples
chunk_association <- function(snps, pos){
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

#' Internal function
#' 
#' Fits a linear model. With data frame and formula.
#' Returns results for variable "Frequency".
#'
#' @param d A data.frame/tibble which must have a "Frequency"
#' column which is the variable of interest. In case of failure
#' or singular value, it returns NA.
#' @param f1 A formula object
#'
#' @return A vector of results.
#'
#' @examples
fit_model <- function(d, f1) {
  #snp <- unique(d$SNP)
  #cat(snp, "\n")
  m1 <- tryCatch(lm(f1, data = d, singular.ok = FALSE),
                 error = function(e){list(coefficients=NA)} )
  if( is.na(m1$coefficients[1]) || is.na(coef(m1)['Frequency']) ){
    r <- rep(NA, 4)
  }else{
    r <- summary(m1)$coefficients["Frequency", ]
  }
  
  names(r) <- c("beta", "SE", "t.value", "p.value")
  return(r)
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
  
  # removing zeroes before log-transformation
  pvector[ pvector == 0 ] <- min(pvector[ pvector > 0 ])
  
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}