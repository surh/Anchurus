# setwd("/cashew/users/sur/exp/fraserv/2021")
library(tidyverse)
library(seqinr)

opts <- commandArgs(trailingOnly = TRUE)
args <- list(snps_info = opts[1],
             snps_alleles = opts[2],
             core_genes = opts[3],
             output = opts[4])
# args <- list(snps_info = "filtered_alleles/MGYG-HGUT-00022/snps_info.txt",
#              snps_alleles = "filtered_alleles/MGYG-HGUT-00022/snps_alleles.txt",
#              core_genes = "filtered_alleles/MGYG-HGUT-00022/core_genes.tsv",
#              output = "core_aln.fasta")


info <- read_tsv(args$snps_info,
                 col_types = cols(site_id = col_character(),
                                  ref_id = col_character(),
                                  ref_pos = col_number(),
                                  gene_id = col_character(),
                                  locus_type = col_character()))
info

# Get list of core genes
core <- read_tsv(args$core_genes,
                 col_types = cols(gene_id = col_character(),
                                  core_gene = col_logical()))
core <- core$gene_id[core$core_gene]

# Get SNP IDs from core genes
snps <- info %>%
  filter(gene_id %in% core) %>%
  select(site_id) %>%
  unlist()

# Get alleles
alleles <- read_tsv(args$snps_alleles,
                    col_types = cols(site_id = col_character(),
                                     .default = col_character()))
alleles <- alleles %>%
  filter(site_id %in% snps)

# Convert to alignment and write
aln <- alleles %>%
  select(-site_id) %>%
  as.list %>%
  map(~ replace_na(.x, "-"))
write.fasta(aln, names = names(aln), file.out = args$output)


