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
    t
  colnames(res) <- site_ids
  return(res)
}

# Original function that used propagate::bigcor to try to attempt
# all by all SNV correlations.
# contig_snv_cors <- function(Dat, contig, depth_thres = 1, snvs = "all"){
#   # contig <- unique(Dat$info$ref_id)[5]
#   # snvs <-"all"
#   # snvs <- "synonymous"
# 
#   cat("Selecting contig", contig, "\n")
#   Dat.contig <- subset_midas(Dat, contig = contig)
# 
#   if(snvs %in% c("synonymous", "non-synonymous")){
#     cat("Selecting", snvs,  "SNVs\n")
#     Dat.contig$info <- HMVAR::determine_snp_effect(Dat.contig$info)
#     selected <- (Dat.contig$info %>%
#       dplyr::filter(!is.na(snp_effect)) %>%
#       dplyr::filter(snp_effect == snvs) %>%
#       dplyr::select(site_id))$site_id
#     Dat.contig <- subset_midas(Dat = Dat.contig,
#                                contig = NULL,
#                                snvs = selected)
#   }else if(snvs == "all"){
#     cat("Using all SNVs\n")
#   }else{
#     stop("ERROR: invalid snvs specification.", call. = TRUE)
#   }
# 
#   freq <- abun2mat(Dat.contig$freq)
#   depth <- abun2mat(Dat.contig$depth)
#   freq[ depth < depth_thres ] <- NA
#   freq <- freq[ ,colSums(is.na(freq)) < nrow(freq) - 1]
# 
#   Cors <- propagate::bigcor(freq, size = min(2000, ncol(freq)),
#                             use = "pairwise.complete.obs")
#   colnames(Cors) <- colnames(freq)
#   row.names(Cors) <- colnames(freq)
# 
#   # Cors[ matrix(c(1:ncol(freq),1:ncol(freq)), ncol = 2) ] <- NA
#   Cors[ which(upper.tri(Cors, diag = TRUE), arr.ind = TRUE) ] <- NA
# 
#   Cors <- Cors[1:ncol(Cors), 1:ncol(Cors)] %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column(var = "site1") %>%
#     tibble::as_tibble() %>%
#     tidyr::pivot_longer(-site1, values_to = "cor") %>%
#     dplyr::filter(!is.na(cor))
# 
#   Pos <- Dat.contig$info %>%
#     select(site_id, ref_pos)
#   Cors <- Cors %>%
#     left_join(Pos, by = c("site1" = "site_id")) %>%
#     select(everything(), pos1 = ref_pos) %>%
#     left_join(Pos, by = c("name" = "site_id")) %>%
#     select(everything(), pos2 = ref_pos) %>%
#     mutate(dist = abs(pos1 - pos2),
#            r2 = cor^2)
# 
#   return(Cors)
# }

# args <- list(midas_dir = "../../../../data/gathered_results/merged.snps/Actinomyces_dentalis_58667/",
#              map = "../../../../data/gathered_results/2019a.hmp.subsite/hmp.subsite_map.txt",
#              depth_thres = 1,
#              min_snvs = 5000,
#              bigcor_dir = "./bigcor",
#              clean = TRUE)
args <- list(midas_dir = opts[1],
            map = opts[2],
            depth_thres = as.numeric(opts[3]),
            min_snvs = as.numeric(opts[4]),
            bigcor_dir = opts[5],
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
options(fftempdir = args$bigcor_dir)
dir.create(args$bigcor_dir)
Cors <- contigs %>%
  map_dfr(~contig_snv_cors(Dat, contig = .,
                           depth_thres = args$depth_thres,
                           snvs = "all"),
          .id = "ref_id")
# Cors

Cors.s <- contigs %>%
  map_dfr(~contig_snv_cors(Dat, contig = .,
                           depth_thres = args$depth_thres,
                           snvs = "synonymous"),
          .id = "ref_id")
# Cors.s

Cors.ns <- contigs %>%
  map_dfr(~contig_snv_cors(Dat, contig = .,
                           depth_thres = args$depth_thres,
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

if(args$clean){
  unlink(args$bigcor_dir, recursive = TRUE)
}

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
