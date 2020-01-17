#!/usr/bin/env Rscript

library(tidyverse)
library(seqinr)
# setwd("/cashew/users/sur/exp/fraserv/2020/today3/gut")

args <- list(gene_coverages = "gene_coverages/",
             filtered_dir = "filtered_alns",
             min_cov = 0.8)

args$gene_coverages <- list.files(args$gene_coverages, full.names = TRUE)
Res <- NULL
for(cov_file in args$gene_coverages){
  # cov_file <- args$gene_coverages[1]
  
  spec <- stringr::str_remove(basename(cov_file),
                              pattern = "\\.gene_coverage\\.txt$")
  cat(spec, "\n")
  filtered_alns <- file.path(args$filtered_dir, paste0(spec, ".concatenated_filtered.aln.fasta"))
  
  # Read coverage
  Cov <- read_tsv(cov_file,
                  col_types = cols(gene = col_character()))
  n_genes <- nrow(Cov)
  n_samples <- ncol(Cov) - 1
  
  # Convert to numeric matrix
  Cov <- Cov %>%
    select(-gene) %>%
    as.matrix
  
  # Filter
  Cov <- Cov[ rowSums(Cov >= args$min_cov) >= (ncol(Cov) * args$min_cov), , drop = FALSE ]
  Cov <- Cov[ , colSums(Cov >= args$min_cov) >= (nrow(Cov) * args$min_cov), drop = FALSE ]
  
  selected_samples <- ncol(Cov)
  selected_genes <- nrow(Cov)
  
  
  if(file.exists(filtered_alns)){
    aln <- seqinr::read.fasta(filtered_alns)
    n_filtered <- length(aln)
    
    n_snvs <- aln %>%
      map_int(length) %>%
      unique
    if(length(n_snvs) != 1){
      stop("ERROR: not all filtered sequences have the same length")
    }
  }else{
    n_filtered <- 0
    n_snvs <- 0
  }
  
  
  Res <- Res %>% 
    bind_rows(tibble(spec = spec,
                     total_genes = n_genes,
                     total_samples = n_samples,
                     selected_genes = selected_genes,
                     selected_samples = selected_samples,
                     filtered_samples = n_filtered,
                     filtered_snvs = n_snvs))
}
Res

write_tsv(Res, "marker_stats.txt")



