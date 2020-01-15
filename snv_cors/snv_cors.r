#!/usr/bin/env Rscript

# (C) Copyright 2019 Sur Herrera Paredes
#
# This file is part of HMVAR.
#
# HMVAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HMVAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HMVAR.  If not, see <http://www.gnu.org/licenses/>.

# setwd("/cashew/users/sur/exp/fraserv/2019/today")
library(tidyverse)
library(HMVAR)
# library(hexbin)
# library(propagate)

#' Subset MIDAS
#'
#' Takes SNV data imported from MIDAS and selects the
#' and returns an object containing a subset of SNVs.
#'
#' At least one of contig and snvs must be specified. In that case,
#' the intersection of both is selected and returned.
#'
#' @param Dat A list with data frame elements info, freq and depth
#' corresponding to the output files of <midas_merge.py snvs>.
#' @param contig A singe contig ID to keep. If NULL, SNVs from all
#' contigs will be kept. If snvs is not NULL, only SNVs both in the
#' contig and in snvs will be kept. Contig name is stored in the
#' ref_id column of Dat$info.
#' @param snvs A character vector of SNV ids to keep. If NULL all
#' SNVs of the given contig will be kept/
#'
#' @return A list with data frame elements info, freq and depth. Similar
#' to the input (Dat)
#' @export
#' @author Sur from Fraser Lab
#'
#' @importFrom magrittr %>%
subset_midas <- function(Dat, contig = NULL, snvs=NULL){
  # contig <- unique(Dat$info$ref_id)[1]

  if(!is.null(snvs) && !is.null(contig)){
    Res <- list(info = Dat$info %>%
                  dplyr::filter(site_id %in% snvs & ref_id == contig))
  }else if(!is.null(snvs) && is.null(contig)){
    Res <- list(info = Dat$info %>%
                  dplyr::filter(site_id %in% snvs))
  }else if(is.null(snvs) && !is.null(contig)){
    Res <- list(info = Dat$info %>%
                  dplyr::filter(ref_id == contig))
  }else{
    stop("ERROR: At least one of snvs and contigs must not be NULL.", call. = TRUE)
  }

  Res$freq <- Dat$freq %>%
    dplyr::filter(site_id %in% Res$info$site_id)
  Res$depth <- Dat$depth %>%
    dplyr::filter(site_id %in% Res$info$site_id)

  return(Res)
}

#' MIDAS abundance to matrix
#'
#' Utility function converts data frame with
#' abundance data into numeric matrix
#'
#' @param tab A data frame with column site_id which will be converted
#' to row  names, and one column per sample.
#'
#' @return A numeric matrix
#' @author Sur from Fraser Lab
abun2mat <- function(tab){
  site_ids <- tab$site_id
  res <- tab %>%
    dplyr::select(-site_id) %>%
    as.matrix
  row.names(res) <- site_ids
  return(res)
}

#' SNV frequency correlations within a contig
#'
#' Calculates pairwise correlation between SNV frequencies
#' for all pairs at a window.
#'
#' @param freqs A n x m numeric matrix of n snvs and m samples with
#' allele frequencies inside.
#' @param positions A n-length vector with the base pair position
#' of each SNV in freqs.
#' @param w_size Non-inclusive maximum distance for a pair of SNVs so
#' that the correlation is calculated. If a SNV is at position x, then
#' all SNVs in the window (x-w_size, x+w_size) will have their correlation
#' calculated.
#' @param circular Eventually handle circular contigs (e.g. bacterial
#' chromosomes, plasmids etc.)
#'
#' @return A data frame with columns pos1, pos2, dist,
#' r2 and n.
#' @author Sur from Fraser Lab
#'
#' @importFrom magrittr %>%
contig_snv_cors <- function(freqs, positions, w_size = 10000, circular = FALSE){

  # Check params
  if(!is.matrix(freqs)){
    stop("ERROR: freqs must be a numeric matrix", call. = TRUE)
  }
  if(length(positions) != nrow(freqs)){
    stop("ERROR: positions must have the same length as the number of rows in freqs", call. = TRUE)
  }

  Res <- NULL
  # Iterate over every
  for(i in 1:nrow(freqs)){
    if((i %% 1000) == 0)
      cat("\tSNV:",i, "\n")

    # Find incremental window
    pos <- positions[i]
    # min_pos <- max(0, pos - w_size)
    max_pos <- pos + w_size

    # Determine indices of snvs to correlate
    # w_snvs <- which(positions != pos & positions > min_pos & positions < max_pos)
    w_snvs <- which(positions > pos & positions < max_pos)

    # Correlate current SNV with other SNVs in window
    rho <- cor(freqs[i,], t(freqs[w_snvs,,drop = FALSE]), use = "pairwise.complete.obs")
    Res <- dplyr::bind_rows(Res,
                            tibble::tibble(pos1 = positions[i],
                                           pos2 = positions[w_snvs],
                                           r2 = rho[1,]^2))
  }

  return(Res)
}

#' SNV frequency correlations for a contig
#'
#' Takes all data from one contig and returns the pairwise
#' correlation between SNVs at a maximum distance.
#'
#' @param Dat A list with data frame elements info, freq and depth
#' corresponding to the output files of <midas_merge.py snvs>.
#' @param contig A singe contig ID to keep. If NULL, SNVs from all
#' contigs will be kept. If snvs is not NULL, only SNVs both in the
#' contig and in snvs will be kept. Contig name is stored in the
#' ref_id column of Dat$info.
#' @param depth_thres Minimum number of reads for a site to be
#' considered covered in a given sample.
#' @param w_size Non-inclusive maximum distance for a pair of SNVs so
#' that the correlation is calculated. If a SNV is at position x, then
#' all SNVs in the window (x-w_size, x+w_size) will have their correlation
#' calculated.
#' @param snvs Either "all", "synonymous" or "non-synonymous".
#'
#' @return A data frame with columns pos1, pos2, dist,
#' r2 and n.
#' @export
#' @author Sur from Fraser Lab
#'
#' @importFrom magrittr %>%
snv_cors <- function(Dat, contig, depth_thres = 1, w_size = 10000, snvs = "all"){

  cat("Selecting contig", contig, "\n")
  Dat.contig <- subset_midas(Dat, contig = contig)

  # Select subset of SNVs
  if(snvs %in% c("synonymous", "non-synonymous")){
    cat("Selecting", snvs,  "SNVs\n")
    Dat.contig$info <- HMVAR::determine_snp_effect(Dat.contig$info)
    selected <- (Dat.contig$info %>%
                   dplyr::filter(!is.na(snp_effect)) %>%
                   dplyr::filter(snp_effect == snvs) %>%
                   dplyr::select(site_id))$site_id
    Dat.contig <- subset_midas(Dat = Dat.contig,
                               contig = NULL,
                               snvs = selected)
  }else if(snvs == "all"){
    cat("Using all SNVs\n")
  }else{
    stop("ERROR: invalid snvs specification.", call. = TRUE)
  }

  # Remove sites and samples not passing depth_thres
  freq <- abun2mat(Dat.contig$freq)
  depth <- abun2mat(Dat.contig$depth)
  freq[ depth < depth_thres ] <- NA
  freq <- freq[ ,colSums(is.na(freq)) < nrow(freq) - 1, drop = FALSE]
  
  if(ncol(freq) > 1){
    # Calculate correlation
    Res <- contig_snv_cors(freqs = freq, positions = Dat.contig$info$ref_pos, w_size = w_size)
  }else{
    Res <- NULL
  }

  return(Res)
}

# args <- list(midas_dir = "../../../../data/gathered_results/merged.snps/Actinomyces_dentalis_58667/",
#              map = "../../../../data/gathered_results/2019a.hmp.subsite/hmp.subsite_map.txt",
#              depth_thres = 1,
#              min_snvs = 5000,
#              clean = TRUE,
#              w_size = 10000)
args <- list(midas_dir = opts[1],
             map = opts[2],
             depth_thres = as.numeric(opts[3]),
             min_snvs = as.numeric(opts[4]),
             w_size = as.numeric(opts[5]),
             clean = as.logical(opts[6]))

# Read data
map <- read_tsv(args$map) %>%
  select(sample = ID, Group)
Dat <- read_midas_data(midas_dir = args$midas_dir, map = map, cds_only = FALSE)

# Select contigs
contigs <- table(Dat$info$ref_id)
contigs <- names(contigs[ contigs >= args$min_snvs ])
contigs <- setNames(contigs, contigs)

# Calculate correlations
Cors <- contigs %>%
  map_dfr(~snv_cors(Dat, contig = .,
                    depth_thres = args$depth_thres,
                    w_size = args$w_size,
                    snvs = "all"),
          .id = "ref_id")
# Cors

Cors.s <- contigs %>%
  map_dfr(~snv_cors(Dat, contig = .,
                    depth_thres = args$depth_thres,
                    w_size = args$w_size,
                    snvs = "synonymous"),
          .id = "ref_id")
# Cors.s

Cors.ns <- contigs %>%
  map_dfr(~snv_cors(Dat, contig = .,
                    depth_thres = args$depth_thres,
                    w_size = args$w_size,
                    snvs = "non-synonymous"),
          .id = "ref_id")
# Cors.ns

# Write output
filename <- "cors_all.txt.gz"
write_tsv(Cors, path = filename)
filename <- "cors_synonymous.txt.gz"
write_tsv(Cors.s, path = filename)
filename <- "cors_nonsynonymous.txt.gz"
write_tsv(Cors.ns, path = filename)

# Plotting will be moved to another script

# filename <- "cors_all_snvs.png"
# png(filename, width = 1000, height = 750, res = 200)
# p1 <- hexbin::hexbin(Cors$dist, Cors$r2)
# hexbin::plot(p1, trans = log10, inv = function(x){10^x}, colramp=terrain.colors)
# dev.off()
#
# filename <- "cors_snynonymous_snvs.png"
# png(filename, width = 1000, height = 750, res = 200)
# p2 <- hexbin::hexbin(Cors.s$dist, Cors.s$r2)
# hexbin::plot(p2, trans = log10, inv = function(x){10^x}, colramp=terrain.colors)
# dev.off()
#
# filename <- "cors_nonsnynonymous_snvs.png"
# png(filename, width = 1000, height = 750, res = 200)
# p3 <- hexbin::hexbin(Cors.ns$dist, Cors.ns$r2)
# hexbin::plot(p3, trans = log10, inv = function(x){10^x}, colramp=terrain.colors)
# dev.off()

# hexbin::hdiffplot(p3,p3)
#
# erodebin1 <- hexbin::erode.hexbin(hexbin::smooth.hexbin(p1))
# erodebin2 <- hexbin::erode.hexbin(hexbin::smooth.hexbin(p3))
# hexbin::hdiffplot(erodebin1, erodebin2)

# x1 <- rnorm(1000)
# x2 <- rnorm(1000)
# y1 <- rnorm(1000)
# y2 <- rnorm(1000)
#
# xbnds <- range(x1,x2)
# ybnds <- range(y1,y2)
#
# p1 <- hexbin::hexbin(x1,y1, xbnds = xbnds, ybnds = ybnds)
# p2 <- hexbin::hexbin(x2,y2, xbnds = xbnds, ybnds = ybnds)
#
# hexbin::hdiffplot(p1, p2)
#
#
# p.diff <- p1
# p.diff@count
#
# p.diff@count - p2@count
#
# length(p.diff@count)
# length(p2@count)
