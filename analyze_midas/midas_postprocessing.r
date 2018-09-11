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
  cat("==== Processing genome ", genome_id, " ====")
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
    return(c(0,0))
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
  return(dim(freqs))
}

#############################

library(plyr)
library(argparser)


## Now need function that goes through midas output directory tree and
## postprocess every genome

meta <- read.table("metadata_qin2012.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(meta)
samples <- meta$ID[ !is.na(meta$Diabetic) ]

missing_value <- NA

input <- "~/micropopgen/exp/2018/2018-09-04.test_merged_snps/merged.snps/"
type <- "dir"
outdir <- "out/"
overwrite <- TRUE
genome_ids_file <- "genome_ids.txt"

if(type == "dir"){
  specdirs <- list.dirs(input, recursive = FALSE, full.names = FALSE)
  genome_ids <- read.table(genome_ids_file, header = TRUE, stringsAsFactors = FALSE)
  
  if(!all(specdirs %in% genome_ids$species)){
    stop("ERROR: species dirs have no ID")
  }
  
  # Create ids lookup table
  ids <- genome_ids$id
  names(ids) <- genome_ids$species
  
  run_table <- data.frame(genome_id = ids[specdirs],
                          freq_file = paste0(input, "/", specdirs, "/snps_freq.txt"),
                          depth_file = paste0(input, "/", specdirs, "/snps_depth.txt"),
                          info_file = paste0(input, "/", specdirs, "/snps_info.txt"),
                          new_freq_file = paste0(outdir, "/", ids[specdirs], ".snps_freq.txt"),
                          new_depth_file = paste0(outdir, "/", ids[specdirs], ".snps_depth.txt"),
                          new_info_file = paste0(outdir, "/", ids[specdirs], ".snps_info.txt"),
                          row.names = ids[specdirs],
                          stringsAsFactors = FALSE)
  # run_table
  
}else if(type == "spec"){
  # A single species dir is passed
  stop("ERROR: type spec not implemented")
}else if(type == "table"){
  # Table of speciesdir, species_name, species_id
  stop("ERROR: type table not implemented")
}

run_table
if(dir.exists(outdir) && !overwrite){
  stop("ERROR: Output dir already exists")
}else if(!dir.exists(outdir)){
  cat("\tCreating directory for output..\n")
  dir.create(outdir)
}
# Call homogenize table
res <- plyr::maply(run_table, homogenize_genome_snps,
                   samples = samples, missing_value = missing_value )
res


