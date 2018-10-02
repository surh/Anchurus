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