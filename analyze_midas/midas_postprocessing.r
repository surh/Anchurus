#' Homogenize genome snps
#' 
#' Function that homogenizes output from MIDAS so that all
#' the samples that are going to be tested are included in the
#' different matrices. It also renames SNPs so that they have
#' a unique name across all genomes
#' 
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
#' 
#' @author Sur Herrera Paredes from Fraser Lab
#' 
#' @export
homogenize_genome_snps <- function(samples, freq_file, depth_file,
                                   info_file, new_freq_file,
                                   new_depth_file, new_info_file,
                                   missing_value = NA){
  
  # Read table of allele frequencies
  freqs <- read.table(freq_file, header = TRUE,
                      sep = "\t", row.names = 1)
  cat("\tProcessing frequencies...\n")
  
  # Select samples that are present
  freqs <- freqs[, intersect(colnames(freqs), samples), drop = FALSE ]
  freqs <- freqs[ rowSums(freqs) > 0, , drop = FALSE ]
  
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
  
  
  return(dim(freqs))
}


library(dplyr)
meta <- read.table("metadata_qin2012.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(meta)

samples <- meta$ID[ !is.na(meta$Diabetic) ]
freq_file <- "snps_freq.txt"
depth_file <- "snps_depth.txt"
info_file <- "snps_info.txt"
genome_id <- "G1"
outdir <- "processed_snps"
missing_value <- NA
new_freq_file <- "processed_snps/G1.snp_freqs.txt"
new_info_file <- "processed_snps/G1.snp_info.txt"
new_depth_file <- "processed_snps/G1.snp_depth.txt"


## Now need function that goes through midas output directory tree and
## postprocess every genome

indir <- "~/micropopgen/exp/2018/2018-09-04.test_merged_snps/merged.snps/"
