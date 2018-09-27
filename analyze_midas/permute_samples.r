#!/usr/bin/env Rscript

library(argparser)
############
#' Process command line arguments
#'
#' @return List of arguments
#' 
#' @author Sur Herrera Paredes
#' 
#' @importFrom argparser arg_parser add_argument parse_args
#' @export
process_arguments <- function(){
  p <- argparser::arg_parser(paste0("Permutes a table of phenotypes the specified ",
                                         "number of times"))
  
  # Positional arguments
  p <- argparser::add_argument(p, "phenotype", help = paste0("File with phenotypes."),
                                    type = "character")
  p <- argparser::add_argument(p, "covariates", help = paste0("File with covariates."),
                               type = "character")
  
  # Optional arguments
  p <- argparser::add_argument(p, "--seed", help = "Seed for random of permutations",
                                    default = 3094229, type = "numeric")
  p <- argparser::add_argument(p, "--nperm", help = "Number of permutations",
                               default = 1000, type = "numeric")
  p <- argparser::add_argument(p, "--prefix", help = "Prefix for output files",
                               default = "phenotype_p", type = "character")
  
  # Read arguments
  args <- argparser::parse_args(parser)
  
  return(args)
}
############

# Define arguments
# args <- list(phenotype = "phenotype.txt",
#              covariates = "covariates.txt",
#              nperm = 10,
#              seed = 3094229,
#              prefix = "phenotype_p")
args <- process_arguments()

pheno <- read.table(args$phenotype, sep = "\t",
                    row.names = 1, header = TRUE)
covs <- read.table(args$covariates, sep = "\t",
                   row.names = 1, header = TRUE)
samples <- colnames(pheno)
if(any(samples != colnames(covs)))
  stop("ERROR: samples don't match")

set.seed(args$seed)
for(i in 1:args$nperm){
  ii <- sample(ncol(pheno), replace = FALSE)
  perm <- pheno[, ii]
  colnames(perm) <- samples
  
  filename <- paste0(args$prefix, i, ".txt")
  write.table(perm, file = filename, col.names = TRUE,
              row.names = TRUE, sep = "\t", quote = FALSE)
  
}

