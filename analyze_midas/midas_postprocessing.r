#!/usr/bin/env Rscript

#' Homogenize genome snps
#' 
#' Function that homogenizes output from MIDAS so that all
#' the samples that are going to be tested are included in the
#' different matrices. It also renames SNPs so that they have
#' a unique name across all genomes
#' 
#' @param genome_id String with genome ID. Used for setting SNP IDs.
#' @param samples Character vector with IDs of samples that will be
#' included in the association testing
#' @param freq_file File with minor allele frequencies. Output from
#' MIDAS.
#' @param depth_file File with sequencing depth per SNP. Output from
#' MIDAS.
#' @param info_file File with information about SNPs. Output from
#' MIDAS.
#' @param new_freq_file Name of new file to be created for minor
#' allele frequencies.
#' @param depth_file Name of new file to be created for sequencing
#' depth per SNP.
#' @param info_file Name of new file to be created for information
#' about SNPs.
#' @param missing_value Value for allele frequencies and snp depths
#' when sample didn't have enough coverage of species.
#' @param new_depth_file 
#' @param new_info_file 
#' 
#' @return The dimensions of the final allele frequency table
#' 
#' @author Sur Herrera Paredes from Fraser Lab
#' 
#' @export
homogenize_genome_snps <- function(genome_id, freq_file, depth_file,
                                   info_file, new_freq_file,
                                   new_depth_file, new_info_file,
                                   samples,
                                   missing_value = NA){
  # i <- 3
  # genome_id <- run_table[i,1]
  # freq_file <- run_table[i,2]
  # depth_file <- run_table[i,3]
  # info_file <- run_table[i,4]
  # new_freq_file <- run_table[i,5]
  # new_depth_file <- run_table[i,6]
  # new_info_file <- run_table[i,7]
  
  # Read table of allele frequencies
  cat("==== Processing genome", genome_id, "====\n")
  freqs <- read.table(freq_file, header = TRUE,
                      sep = "\t", row.names = 1)
  cat("\tProcessing frequencies...\n")
  
  # Select samples that are present
  freqs <- freqs[, intersect(colnames(freqs), samples), drop = FALSE ]
  freqs <- freqs[ rowSums(freqs) > 0, , drop = FALSE ]
  
  # Check if there are samples and SNPs to process
  if(any(dim(freqs) == 0)){
    cat("\tNo samples or SNPs found for this genome.\n")
    cat("==== DONE ====\n")
    return(data.frame(nsamples = ncol(freqs),
                      nsnes = nrow(freqs)))
  }
  
  # Obtain sequencing depths and filter according to allele freqs
  depth <- read.table(depth_file, header = TRUE,
                      sep = "\t", row.names = 1)
  cat("\tProcessing dephts...\n")
  depth <- depth[ row.names(freqs), colnames(freqs), drop = FALSE ]
  
  # Read and filter SNP info table
  info <- read.table(info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat("\tProcessing info...\n")
  row.names(info) <- as.character(info$site_id)
  info <- info[ row.names(freqs), , drop = FALSE]
  
  # Change names of SNPs
  cat("\tRenaming...\n")
  new_ids <- paste(genome_id, row.names(freqs), sep = ".")
  row.names(freqs) <- new_ids
  row.names(depth) <- new_ids
  info$site_id <- new_ids
  
  # Add samples that were missing because species not above threshold
  cat("\tAdding missing samples...\n")
  missing_samples <- matrix(missing_value, nrow = nrow(freqs), ncol = length(samples) - ncol(freqs),
                            dimnames = list(row.names(freqs),
                                            setdiff(samples, colnames(freqs))))
  freqs <- cbind(freqs, missing_samples)
  freqs <- freqs[ , samples ]
  depth <- cbind(depth, missing_samples)
  depth <- depth[ , samples ]
  
  # Write output
  cat("\tWriting output...\n")
  write.table(freqs, new_freq_file, sep = "\t",
              quote = FALSE, col.names = NA, row.names = TRUE)
  write.table(depth, new_depth_file, sep = "\t",
              quote = FALSE, col.names = NA, row.names = TRUE)
  write.table(info, new_info_file, sep = "\t",
              quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  cat("==== DONE ====\n")
  return(data.frame(nsamples = ncol(freqs),
                    nsnes = nrow(freqs)))
}

#' Read arguments
#' 
#' Internal function
#'
#' @return
#' @examples
#' 
#' @importFrom argparser arg_parser add_argument parse_args
process_arguments <- function(){
  parser <- argparser::arg_parser(paste0("Post-process MIDAS mege SNPs output. ",
                                         "Create tables with homogeneous, samples ",
                                         "and descriptions"))
  
  # Positional arguments
  parser <- argparser::add_argument(parser, "type", help = paste0("Type of input to be passed. ",
                                                                  "Must be dir, specdir or table."),
                                    type = "character")
  parser <- argparser::add_argument(parser, "input", help = paste0("Input dir or file. Processed ",
                                                                   "according to 'type'"),
                                    type = "character")
  parser <- argparser::add_argument(parser, "samples", help = "File containing the samples to keep",
                                    type = "character")
  
  # Optional arguments
  parser <- argparser::add_argument(parser, "--outdir", help = "Output directory",
                                    default = "out/", type = "character")
  parser <- argparser::add_argument(parser,"--missing_value", help = "How to encode missing value",
                                    default = NA, type = "numeric")
  parser <- argparser::add_argument(parser, "--overwrite", help = "Whether to overwrite exisiting files",
                                    flag = TRUE, default = FALSE)
  parser <- argparser::add_argument(parser, "--genome_ids", help = paste0("For type dir this must be a file ",
                                                                          "mapping genome names to unique IDs.\n",
                                                                          "For type specdir thus must be a single ",
                                                                          "string corresponding to the genome id."),
                                    default = NULL, type = "character")
  parser <- argparser::add_argument(parser, "--results", help = "File to write post-processing results.",
                                    default = "postprocessing_results.txt", type = "character")
  
  # Read arguments
  args <- argparser::parse_args(parser)
  
  return(args)
}
#############################

library(plyr)
library(argparser)

# Read arguments
args <- process_arguments()

# Get samples
samples <- read.table(args$samples, stringsAsFactors = FALSE)$V1

# Create run table
if(args$type == "dir"){
  if(is.null(args$genome_ids))
    stop("ERROR: genome_ids cannot be NULL for type dir")
  
  specdirs <- list.dirs(args$input, recursive = FALSE, full.names = FALSE)
  genome_ids <- read.table(args$genome_ids, header = TRUE, stringsAsFactors = FALSE)
  
  if(!all(specdirs %in% genome_ids$species)){
    stop("ERROR: species dirs have no ID")
  }
  
  # Create ids lookup table
  ids <- genome_ids$id
  names(ids) <- genome_ids$species
  
  run_table <- data.frame(genome_id = ids[specdirs],
                          freq_file = paste0(args$input, "/", specdirs, "/snps_freq.txt"),
                          depth_file = paste0(args$input, "/", specdirs, "/snps_depth.txt"),
                          info_file = paste0(args$input, "/", specdirs, "/snps_info.txt"),
                          new_freq_file = paste0(args$outdir, "/", ids[specdirs], ".snps_freq.txt"),
                          new_depth_file = paste0(args$outdir, "/", ids[specdirs], ".snps_depth.txt"),
                          new_info_file = paste0(args$outdir, "/", ids[specdirs], ".snps_info.txt"),
                          row.names = ids[specdirs],
                          stringsAsFactors = FALSE)
}else if(args$type == "specdir"){
  # A single species dir is passed
  # stop("ERROR: type specdir not implemented")
  if(length(args$genome_ids) != 1)
    stop("ERROR: genome_ids must be a single character string with the specdir type.")
  
  # Create run_table
  run_table <- data.frame(genome_id = args$genome_ids,
                          freq_file = paste0(args$input, "/snps_freq.txt"),
                          depth_file = paste0(args$input, "/snps_depth.txt"),
                          info_file = paste0(args$input, "/snps_info.txt"),
                          new_freq_file = paste0(args$outdir, "/snps_freq.txt"),
                          new_depth_file = paste0(args$outdir, "/snps_depth.txt"),
                          new_info_file = paste0(args$outdir, "/snps_info.txt"),
                          row.names = args$genome_ids,
                          stringsAsFactors = FALSE)
  
}else if(args$type == "table"){
  # Table of speciesdir, species_name, species_id
  stop("ERROR: type table not implemented")
}else{
  stop("ERROR: not recognized type. Must be dir, specdir or table")
}

# Check of output dir exists and create if required
if(dir.exists(args$outdir) && !args$overwrite){
  stop("ERROR: Output dir already exists")
}else if(!dir.exists(args$outdir)){
  cat("\tCreating directory for output..\n")
  dir.create(args$outdir)
}

# Call homogenize table
res <- plyr::mdply(run_table, homogenize_genome_snps,
                   samples = samples, missing_value = args$missing_value )

# Write results
if(!is.null(args$results))
  write.table(res, args$results, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
