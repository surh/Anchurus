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

library(argparser)

#################### FUNCTIONS ##################
process_arguments <- function(){
  p <- arg_parser(paste("produce alignments from metagenomes"))
  
  # Positional arguments
  p <- add_argument(p, "midas_dir",
                    help = paste("Directory with merged MIDAS SNVs."),
                    type = "character")
  
  p <- add_argument(p, "genomes_dir",
                    type = "character",
                    help = paste("Path to directory containing genome directories.",
                                 "Must have one subdirectory per species, and the",
                                 "species name must match the midas_dir species name."))
  
  
  # Optional arguments
  p <- add_argument(p, "--min_cov",
                     help = paste("Minumum coverage required. Must be a proportion",
                                  "between 0 and 1. For any given samples only genes",
                                  "that have at least <min_cov> of their snps sequenced",
                                  "are considered. Further, only genes that are detected",
                                  "at that coverage in at least <min_cov> of the samples",
                                  "are considered. Finally, only samples when at least",
                                  "<min_cov> of the genes are covered above the threshold",
                                  "are kept."),
                     type = "numeric",
                     default = 0.8)
  p <- add_argument(p, "--map_file",
                    help = paste("Name of file with map. Must have ID and Group columns"),
                    default = "map.txt",
                    type = "character")
  p <- add_argument(p, "--outdir",
                    help = paste("Directory for output"),
                    default = "output/",
                    type = "character")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(args$min_cov < 0 || args$min_cov > 1){
    stop("ERROR: --min_cov must be in [0,1]", call. = TRUE)
  }
  # Adding parameters that user cannot modify.
  args$depth_thres <- 1
  args$freq_thres <- 1
  args$keep_last_codon <- TRUE
  args$missing_as <- "gap"
  
  return(args)
}

# Read data
read_data <- function(spec, midas_dir, genomes_dir){
  genome_fasta <- seqinr::read.fasta(file.path(genomes_dir, spec, 'genome.fna.gz'))
  genome_feats <- readr::read_tsv(file.path(genomes_dir, spec, 'genome.features.gz'),
                                  col_types = readr::cols(.default = readr::col_character(),
                                                          start = readr::col_integer(),
                                                          end = readr::col_integer()))
  
  # Select CDS genes
  genome_feats <- genome_feats %>%
    filter(gene_type == "CDS")
  
  # Read midas snv data
  Dat <- read_midas_data(midas_dir = midas_dir,
                         map = map,
                         genes = genome_feats$gene_id)
  
  
  return(list(genome_fasta = genome_fasta,
              genome_feats = genome_feats,
              midas = Dat))
}
##################################

args <- process_arguments()
# args <- list(depth_thres = 1,
#              freq_thres =  0.5,
#              min_cov =  0.8,
#              keep_last_codon = TRUE,
#              outdir = 'output/',
#              genomes_dir = "/home/sur/micropopgen/data/genomes/midas_db_v1.2/hmp.subsite/",
#              map_file = "midas/map.txt",
#              midas_dir = "midas/merged.snps/Veillonella_parvula_57794/",
#              missing_as = 'gap')

# Load rest of libraries
library(tidyverse)
library(HMVAR)
library(seqinr)

# Create ouput dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

# Read map
map <- read_tsv(args$map_file,
                col_types = cols(.default = col_character())) %>%
  select(sample = ID, Group)

# Process midas_dir
# Eventually a directory of midas dirs can be passed.
for(midas_dir in args$midas_dir){
  # midas_dir <- args$midas_dir[1]
  spec <- basename(midas_dir)
  cat(spec, "\n")
  
  # Read data
  cat("\tReading data...\n")
  Dat <- read_data(spec = spec, midas_dir = midas_dir, genomes_dir = args$genomes_dir)
  # map %>% filter(sample %in% colnames(Dat$midas$freq)) %>% select(Group) %>% table
  
  if(ncol(Dat$midas$freq) > 2){
    # At least two samples for analysis
    
    # For each gene, extract MSA
    cat("\nExtracting MSAs...\n")
    Res <- Dat$midas$info %>%
      # head(30000) %>%  ##!! JUST FOR TESTING
      # head(1000) %>%  ##!! JUST FOR TESTING
      split(.$gene_id) %>%
      map(function(i, freq, depth, map, depth_thres, freq_thres, min_cov,
                   genome_feats, missing_as, genome_fasta){
        
        # Get gene information
        genome_feats <- genome_feats %>%
          filter(gene_id == unique(i$gene_id))
        gene <- genome_feats$gene_id
        gene_start <- genome_feats$start
        gene_end <- genome_feats$end
        gene_ref_id <- genome_feats$scaffold_id
        gene_strand <- genome_feats$strand
        
        # Get SNVs from gene
        f <- freq %>% dplyr::filter(site_id %in% i$site_id)
        d <- depth %>% dplyr::filter(site_id %in% i$site_id)
        
        # Calculate coverage and discard samples with uncovered gene
        # gene_cov if the proportion of SNVs in a gene covered per sample
        gene_cov <- colSums((d %>% select(-site_id)) >= depth_thres) / nrow(i)
        if(sum(gene_cov >= min_cov) < min_cov*length(gene_cov)){
          # If less than min_cov of genes have gene_cog >= min_cov, then discard that gene.
          gene_cov <- tibble(gene = gene, sample = names(gene_cov), coverage = gene_cov)
          return(list(Coverage = gene_cov, aln = NULL))
        }
        
        # Select samples that have at least min_cov coverage of gene
        d <- d %>% select(site_id, names(which(gene_cov >= min_cov)))
        f <- f %>% select(site_id, names(which(gene_cov >= min_cov)))
        gene_cov <- tibble(gene = gene, sample = names(gene_cov), coverage = gene_cov)
        
        # Match freq and depth
        d <- d %>% tidyr::gather(key = "sample", value = "depth",
                                 -site_id)
        f <- f %>% tidyr::gather(key = "sample", value = "freq",
                                 -site_id)
        dat <- d %>%
          dplyr::left_join(f, by = c("site_id",
                                     "sample")) %>%
          dplyr::filter(depth >= depth_thres)
        dat <- dat %>% dplyr::left_join(map, by = "sample")
        dat <- dat %>% dplyr::left_join(i, by = "site_id")
        
        # Extract sequence fragment
        gene_seq <- HMVAR:::get_sequence_fragment(genome_fasta,
                                                  ref_id = gene_ref_id,
                                                  start = gene_start,
                                                  end = gene_end)
        
        dat <- dat %>%
          mutate(ref_pos = ref_pos - gene_start + 1) %>%
          HMVAR:::assign_allele(depth_thres = depth_thres, freq_thres = freq_thres, sequence = TRUE, na_rm = FALSE)
        
        aln <- gene_snv_aln(snps = dat, seq = gene_seq,
                            missing_as = missing_as,
                            strand = gene_strand,
                            keep_last_codon = TRUE)
        
        return(list(Coverage = gene_cov, aln=aln$seq))
      }, freq=Dat$midas$freq, depth=Dat$midas$depth, map=map,
      depth_thres=args$depth_thres, freq_thres=args$freq_thres, min_cov=args$min_cov,
      genome_feats=Dat$genome_feats, missing_as=args$missing_as,
      genome_fasta=Dat$genome_fasta)
    
    # Produce and write coverage table
    cat("\tProducing gene coverage table...\n")
    cov <- Res %>%
      map_dfr(~ .$Coverage) %>%
      spread(sample, coverage)
    filename <- file.path(args$outdir, paste0(spec, ".gene_coverage.txt"))
    readr::write_tsv(cov, path = filename)
    
    # Select genes covered at least min_cov in at least min_cov of the samples
    cat("\tSelecting genes...\n")
    genes <- cov$gene
    cov <- cov %>%
      select(-gene) %>%
      as.matrix
    row.names(cov) <- genes
    cov <- cov[ rowSums(cov >= args$min_cov) >= (ncol(cov) * args$min_cov), ]
    cov <- cov[ , colSums(cov >= args$min_cov) >= (nrow(cov) * args$min_cov) ]
    
    # Write selected alignments
    cat("\tWriting alns...\n")
    n_seqs <- row.names(cov) %>%
      purrr::map(function(gene, Res, samples, outdir){
        aln <- Res[[gene]]$aln
        aln <- aln[intersect(names(aln), samples)]
        filename <- file.path(outdir, paste0(gene, ".aln.fasta"))
        seqinr::write.fasta(sequences = aln, names = names(aln), file.out = filename)
        return(length(aln))
      }, Res = Res, samples = colnames(cov), outdir = args$outdir)
    
  }else{
    cat("No samples for", spec, "\n")
    Res <- NULL
  }
}

