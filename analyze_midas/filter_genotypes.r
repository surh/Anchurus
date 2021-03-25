#!/usr/bin/env Rscript
library(tidyverse)
library(HMVAR)
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Filters samples and genotype calls"))

  # Positional arguments
  p <- add_argument(p, "genotypes",
                    help = paste("A genotype file in matrix format",
                                 "Must be a tab delimited file with",
                                 "a column 'site_id', and one additional",
                                 "column per sample. A dot is interpreted",
                                 "as missing data."),
                    type = "character")
  p <- add_argument(p, "info",
                    help = paste("A file with SNV information. Must match the",
                                 "format of snps_info.txt created by MIDAS."),
                    type = "character")
  p <- add_argument(p, "map",
                    help = paste("Mapping file. Must have columns 'sample',",
                                 "and group."),
                    type = "character")

  # Optional arguments
  p <- add_argument(p, "--outdir",
                    help = paste("Directory to write output"),
                    type = "character",
                    default = "output/")
  p <- add_argument(p, "--min_samples_per_group",
                    help = paste("Expects samples to belong to one of two",
                                 "groups. Minimun number of samples per",
                                 "group to keep whole species and SNVs.",
                                 "If fewer samples of either group pass the",
                                 "threshold there won't be a final genotypes",
                                 "file."),
                    type = "integer",
                    default = 5)
  p <- add_argument(p, "--min_snv_prop_per_sample",
                    help = paste("Minimum proportion of total SNVs in a",
                                 "sample."),
                    type = "numeric",
                    default = 0.5)
  p <- add_argument(p, "--min_core_gene_prev",
                    help = paste("Minimum proportion of samples (prevalence)",
                                 "where a gene is present in order to be",
                                 "'core'"),
                    type = "numeric",
                    default = 0.8)
  p <- add_argument(p, "--min_core_gene_cov",
                    help = paste("Minimum fraction of a gene's SNVs that are",
                                 "covered for that gene to be present."),
                    type = "numeric",
                    default = 0.8)
  p <- add_argument(p, "--min_core_genes",
                    help = paste("Minimum proportion of core genes in a",
                                 "sample."),
                    type = "numeric",
                    default = 0.8)
  p <- add_argument(p, "--process_info",
                    help = paste("Indicates whether the info file should also",
                                 "be proessed. Must be 'yes' or 'no'"),
                    type = "character",
                    default = "yes")

  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)

  # Process arguments
  if(args$min_snv_prop_per_sample < 0 | args$min_snv_prop_per_sample < 0)
    stop("ERROR: min_snv_prop_per_sample must be in [0,1]", call. = TRUE)
  if(args$min_core_gene_prev < 0 | args$min_core_gene_prev < 0)
    stop("ERROR: min_core_gene_prev must be in [0,1]", call. = TRUE)
  if(args$min_core_gene_cov < 0 | args$min_core_gene_cov < 0)
    stop("ERROR: min_core_gene_cov must be in [0,1]", call. = TRUE)
  if(args$min_core_genes < 0 | args$min_core_genes < 0)
    stop("ERROR: min_core_genes must be in [0,1]", call. = TRUE)
  if(!(args$process_info %in% c("yes", "no")))
    stop("ERROR: process_info must be 'yes' or 'no'", call. = TRUE)

  args$process_info <-  purrr::set_names(x = c(TRUE, FALSE),
                                         nm = c("yes", "no"))[args$process_info]

  return(args)
}

check_groups <- function(vec, ngroups = 2, minsize = 5){
  tab <- table(vec[!is.na(vec)])
  length(tab) == ngroups && all(tab >= minsize)
}

args <- process_arguments()
# args <- list(genotypes = "/home/sur/micropopgen/exp/2020/today/snps_alleles.txt",
#              info = "/home/sur/micropopgen/exp/2020/today/MGYG-HGUT-04165/snps_info.txt",
#              map = "/home/sur/micropopgen/exp/2020/today/example_map.txt",
#              outdir = "/home/sur/micropopgen/exp/2020/today/output",
#              min_samples_per_group = 2,
#              min_snv_prop_per_sample = 0.5,
#              min_core_gene_prev = 0.8,
#              min_core_gene_cov = 0.8,
#              min_core_genes = 0.8,
#              process_info = TRUE)


cat("Reading data...\n")
# Read file with grouping variable
map <- read_tsv(args$map,
                col_types = cols(sample = col_character(),
                                 group = col_character()))

# Read info file from midas merge
info <- readr::read_tsv(args$info,
                        col_types = readr::cols(.default = readr::col_character(),
                                                ref_pos = readr::col_number(), count_samples = readr::col_number(),
                                                count_a = readr::col_number(), count_c = readr::col_number(),
                                                count_g = readr::col_number(), count_t = readr::col_number()),
                        na = "NA")

# Read genotype calls in matrix format
geno <- read_tsv(file = args$genotypes,
                 col_types = cols(site_id = col_character(),
                                  .default = col_character()),
                 na = c(".", "NA", ""))

# Select samples in map
geno <- geno[ ,colnames(geno) %in% c("site_id", map$sample) ]
map <- map %>%
  filter(sample %in% colnames(geno)[-1])
if(!check_groups(map$group, ngroups = 2,
                 minsize = args$min_samples_per_group)){
  cat("Not enough initial samples with enough SNVs...")
  quit(save = "no")
}

# Filter samples
cat("Filtering samples without enough SNVs\n")
n_snvs <- nrow(geno)
# n_samples <- ncol(geno) - 1
n_missing_per_sample <- colSums(is.na(geno[, -1]))
# n_missing_per_snv <- rowSums(is.na(geno[, -1]))
# 1 - n_missing_per_sample / nrow(geno)
# ncol(geno) - 1 - n_missing_per_snv
selected_samples <- names(n_missing_per_sample)[1 - n_missing_per_sample / n_snvs > args$min_snv_prop_per_sample]
geno <- geno %>%
  dplyr::select(site_id, all_of(selected_samples))
map <- map %>%
  filter(sample %in% colnames(geno)[-1])
if(!check_groups(map$group, ngroups = 2,
                 minsize = args$min_samples_per_group)){
  cat("Not enough samples per group with enough SNVs...")
  quit(save = "no")
}

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}else{
  stop("ERROR: output directory already exists", call. = TRUE)
}

# Get SNV coverage fraction per gene
cat("Calculating fraction of SNV coverage per gene...\n")
gene_cov <- info %>%
  dplyr::select(site_id, locus_type, gene_id) %>%
  dplyr::filter(locus_type == "CDS") %>%
  dplyr::select(-locus_type) %>%
  dplyr::inner_join(geno, by = "site_id") %>%
  dplyr::select(-site_id) %>%
  # filter(gene_id == "GUT_GENOME258268_00001") %>%
  # filter(gene_id == "GUT_GENOME258268_00002") %>%
  # select(-gene_id) %>% apply(2, function(x){sum(!is.na(x)) / length(x)})
  dplyr::group_by(gene_id) %>%
  summarise(across(.fns = function(x){sum(!is.na(x)) / length(x)}),
            .groups = 'drop')
write_tsv(gene_cov, file.path(args$outdir, "gene_coverage.tsv"))

# Find core genes
cat("Finding core genes...\n")
core_genes <- tibble(gene_id = gene_cov$gene_id,
                     core_gene = apply(gene_cov[,-1], 1, function(x, min_core_gene_prev = 0.8,
                                                                  min_core_gene_cov = 0.8){
                       (sum(x >= min_core_gene_cov) / length(x)) >= min_core_gene_prev
                     }, min_core_gene_prev = args$min_core_gene_prev,
                     min_core_gene_cov = args$min_core_gene_cov))
# core_genes
# table(core_genes$core_gene)
write_tsv(core_genes, file.path(args$outdir, "core_genes.tsv"))

# Find genotypes with enough core genes
cat("Calculating core genes per genome...\n")
# gene_cov[core_genes$core_gene,]
core_genes_per_genome <- gene_cov[core_genes$core_gene,] %>%
  select(-gene_id) %>%
  # group_by(gene_id) %>%
  summarise(across(.fns = function(x, min_core_gene_cov = 0.8){
    sum(x >= min_core_gene_cov) / length(x)
  }, min_core_gene_cov = args$min_core_gene_prev))
# core_genes_per_genome
core_genes_per_genome <- purrr::set_names(x = as.numeric(core_genes_per_genome),
                                          nm = colnames(core_genes_per_genome))
# core_genes_per_genome
write_tsv(tibble::tibble(sample = names(core_genes_per_genome),
                         core_genes_per_genome = as.numeric(core_genes_per_genome)),
          file.path(args$outdir, "core_genes_per_sample.tsv"))
geno <- geno[ ,c(TRUE, core_genes_per_genome >= args$min_core_genes)]
map <- map %>%
  filter(sample %in% colnames(geno)[-1])
if(!check_groups(map$group, ngroups = 2,
                 minsize = args$min_samples_per_group)){
  cat("Not enough samples per group with enough core genes...\n")
  quit(save = "no")
}

# Filter constant SNVs
cat("Removing constant positions...\n")
bi_snvs <- apply(geno[, -1], 1, function(g){
  length(unique(g[!is.na(g)]))
}) == 2
geno <- geno[bi_snvs, ]
# geno

# Filter SNVs that are not in enough samples
cat("Filtering out SNVs without enough samples per group...\n")
snvs_to_keep <- apply(geno[,-1], 1, function(vec, meta, min_samples_per_group = 5){
  vec <- vec[!is.na(vec)]
  tab <- table(meta[ names(vec) ])

  length(tab) == 2 && all(tab >= min_samples_per_group)
}, meta = purrr::set_names(x = map$group, nm = map$sample),
min_samples_per_group = args$min_samples_per_group)
# sum(snvs_to_keep)
geno <- geno[snvs_to_keep, ]
# geno
cat("Writing final set of genotypes...\n")
write_tsv(geno, file.path(args$outdir, "snps_alleles.txt"))

if(args$process_info){
  # info from selected snvs
  cat("Processing info file...\n")
  info <- geno %>%
    select(site_id) %>%
    left_join(info, by = 'site_id') %>%
    select(site_id, ref_id,ref_pos, ref_allele, major_allele,
           minor_allele,
           locus_type, gene_id, snp_type, site_type, amino_acids) %>%
    # select(locus_type) %>% table
    HMVAR::determine_snp_effect() %>%
    HMVAR::determine_substitution_type()
  # info
  cat("Wrting processed info...\n")
  write_tsv(info, file.path(args$outdir, "snps_info.txt"))
}
