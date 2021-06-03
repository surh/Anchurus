#!/usr/bin/env Rscript

# (C) Copyright 2021 Sur Herrera Paredes
# This file is part of Anchurus.
# 
# Anchurusis free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Anchurusis distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Anchurus.  If not, see <https://www.gnu.org/licenses/>.

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Plot heatmap of gene coverage per sample for QP'd",
                        "samples."))
  
  # Positional arguments
  p <- add_argument(p, "input",
                    help = paste("Directory with filtered SNPs."),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--type",
                     help = paste("Either 'multi' or 'single'"),
                     type = "character",
                     default = "single")
  p <- add_argument(p, "--output",
                    help = paste("Filename for single output."),
                    type = "character",
                    default = "heatmap.jpeg")
  p <- add_argument(p, "--outdir",
                    help = paste("Output directory for multi"),
                    type = "character",
                    default = "output")
  p <- add_argument(p, "--taxonomy",
                    help = "File with genome taxonomy",
                    type = "character",
                    default = NULL)
  p <- add_argument(p, "--width",
                    help = paste("Plot width in inches"),
                    type = "numeric",
                    default = 6)
  p <- add_argument(p, "--height",
                    help = paste("Plot width in inches"),
                    type = "numeric",
                    default = 4)

                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(!(args$type %in% c("single", "multi"))){
    stop("ERROR: --type must be 'single' or 'multi'", call. = TRUE)
  }
  
  if(args$type == "multi"){
    args$input <- list.dirs(args$input, full.names = TRUE, recursive = FALSE)
    if(length(args$input) == 0){
      stop("ERROR: There are no subdirectories in input", call. = TRUE)
    }
    
  }
    
  return(args)
}

plot_heatmap <- function(gene_coverage_file,
                         core_genes_file,
                         filtered_file,
                         # outfile = "heatmap.jpeg",
                         width = 6,
                         height = 4,
                         taxonomy = NULL){
  gene_cov <- read_tsv(gene_coverage_file,
                       col_types = cols(gene_id = col_character(),
                                        .default = col_number()))
  core_genes <- read_tsv(core_genes_file,
                         col_types = cols(gene_id = col_character(),
                                          core_gene = col_logical()))
  final_samples <- read_lines(filtered_file, n_max = 1) %>%
    str_split(pattern = "\t")
  final_samples <- final_samples[[1]][-1]
  
  genome_id <- core_genes$gene_id %>%
    str_split(pattern = "_") %>%
    map_chr(function(x){
      paste(x[1], x[2], sep = "_")
    }) %>%
    unique
  # print(genome_id)
  
  # Include taxonomy if present
  if(!is.null(taxonomy)){
    taxonomy <- taxonomy$Lineage[ taxonomy$Genome == genome_id]
    taxonomy <- str_split(taxonomy, ";")[[1]]
    # print(taxonomy)
    for(i in rev(1:7)){
      if(!(taxonomy[i] %in% c("s__","g__","f__","o__","c__", "p__", "d__"))){
        genome_tax <- taxonomy[i]
        break
      }
    }
  }else{
    genome_tax <- NULL
  }
  
  gene_cov %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "frac_covered") %>%
    mutate(core = "accessory",
           selected = "discarded") %>%
    mutate(core = replace(core, gene_id %in% core_genes$gene_id[core_genes$core_gene], "core"),
           selected = replace(selected, sample %in% final_samples, "selected")) %>%
    ggplot(aes(x = gene_id, y = sample)) +
    facet_grid(selected ~ core, scales = "free", space = "free") +
    geom_tile(aes(fill = frac_covered)) +
    scale_fill_gradient(low = "#fff5f0", high = "#67000d",
                        name = "Gene coverage\nbreadth") +
    ggtitle(label = genome_tax) +
    theme(axis.text = element_blank(),
          panel.background = element_blank(),
          legend.position = "top",
          legend.background = element_blank(),
          legend.key = element_blank())
  
}


# setwd("/home/sur/micropopgen/exp/2021/today")
# args <- list(gene_coverage = "filtered_alleles/gene_coverage.tsv",
#              core_genes = "filtered_alleles/core_genes.tsv",
#              filtered_alleles = "filtered_alleles/snps_alleles.txt",
#              tyoe = "single",
#              taxonomy = "genomes-all_lineages.tsv.gz",
#              output = "heatmap.jpeg",
#              outdir = "output",
#              width = 6,
#              height = 4)
args <- process_arguments()

library(tidyverse)
library(HMVAR)

if(args$type == "multi"){
  # Prepare output dir
  if(!dir.exists(args$outdir)){
    dir.create(args$outdir)
  }
}

if(!is.na(args$taxonomy)){
  taxonomy <- read_tsv(args$taxonomy,
                       col_types = cols(Genome = col_character(),
                                        Lineage = col_character()))
}else{
  taxonomy <- NULL
}

for(d in args$input){
  
  if(args$type == "multi"){
    spec <- basename(d)
    args$output = file.path(args$outdir, paste0(spec, ".jpeg"))
  }
  filtered_file <- file.path(d, "snps_alleles.txt")
  if(!file.exists(filtered_file)){
    warning(paste("No snps_alleles.txt file for directory", spec, "skipping."), immediate. = TRUE)
    next
  }
  p1 <- plot_heatmap(gene_coverage_file = file.path(d, "gene_coverage.tsv"),
                     core_genes_file = file.path(d, "core_genes.tsv"),
                     filtered_file = file.path(d, "snps_alleles.txt"),
                     width = args$width,
                     height = args$height,
                     taxonomy = taxonomy)
  ggsave(args$output, p1, width = args$width, height = args$height, dpi = 150)
}


# 
# tax <- read_tsv(args$taxonomy)
# tax <- tax$Lineage[ tax$Genome == genome_id]
# tax <- str_split(tax, ";")[[1]]
# for(i in rev(1:7)){
#   if(!(tax[i] %in% c("s__","g__","f__","o__","c__", "p__", "d__"))){
#     genome_tax <- tax[i]
#     break
#   }
# }
# genome_tax
# 
# # p1
# filename <- args$output
# ggsave(filename, p1, width = args$width, height = args$height, dpi = 150)




