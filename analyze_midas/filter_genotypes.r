library(tidyverse)
library(HMVAR)

args <- list(genotypes = "/home/sur/micropopgen/exp/2020/today/snps_alleles.txt",
             info = "/home/sur/micropopgen/exp/2020/today/MGYG-HGUT-04165/snps_info.txt",
             map = "/home/sur/micropopgen/exp/2020/today/example_map.txt",
             outdir = "output",
             min_samples_per_group = 2,
             min_snv_prop_per_sample = 0.5,
             min_core_gene_prev = 0.8,
             min_core_gene_cov = 0.8,
             min_core_genes = 0.8)



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
s
# Select samples in map
geno <- geno[ ,colnames(geno) %in% c("site_id", map$sample) ]
map <- map %>%
  filter(sample %in% colnames(geno)[-1])
map

# Filter samples
n_snvs <- nrow(geno)
n_samples <- ncol(geno) - 1

n_missing_per_sample <- colSums(is.na(geno[, -1]))
n_missing_per_snv <- rowSums(is.na(geno[, -1]))

1 - n_missing_per_sample / nrow(geno)
ncol(geno) - 1 - n_missing_per_snv

selected_samples <- names(n_missing)[1 - n_missing_per_sample / n_snvs > args$min_snv_prop_per_sample]
geno <- geno %>%
  dplyr::select(site_id, selected_samples)




# Here I would need to determine core genes and 
# samples with corresponding strains

# Get SNV coverage fraction per gene
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
gene_cov

# Find core genes
core_genes <- tibble(gene_id = gene_cov$gene_id,
       core_gene = apply(gene_cov[,-1], 1, function(x, min_core_gene_prev = 0.8,
                                                    min_core_gene_cov = 0.8){
         (sum(x >= min_core_gene_cov) / length(x)) >= min_core_gene_prev 
       }, min_core_gene_prev = args$min_core_gene_prev,
       min_core_gene_cov = args$min_core_gene_cov))
core_genes
table(core_genes$core_gene)

# Find genotypes with enough core genes
gene_cov[core_genes$core_gene,]
core_genes_per_genome <- gene_cov[core_genes$core_gene,] %>%
  select(-gene_id) %>%
  # group_by(gene_id) %>%
  summarise(across(.fns = function(x, min_core_genes_cov = 0.8){
    sum(x >= min_core_genes_cov) / length(x)
  }))

core_genes_per_genome
core_genes_per_genome <- purrr::set_names(x = as.numeric(core_genes_per_genome),
                                          nm = colnames(core_genes_per_genome))
core_genes_per_genome
geno <- geno[ ,c(TRUE, core_genes_per_genome >= args$min_core_genes)]
geno

# Filter constant SNVs
bi_snvs <- apply(geno[, -1], 1, function(g){
  length(unique(g[!is.na(g)]))
}) == 2
geno <- geno[bi_snvs, ]
geno


# Filter SNVs that are not in enough samples
snvs_to_keep <- apply(geno[,-1], 1, function(vec, meta, min_samples_per_group = 5){
  vec <- vec[!is.na(vec)]
  tab <- table(meta[ names(vec) ])

  length(tab) == 2 && all(tab >= min_samples_per_group)
  }, meta = meta,
  min_samples_per_group = args$min_samples_per_group)
sum(snvs_to_keep)
geno <- geno[snvs_to_keep, ]
geno

# info from selected snvs
info <- geno %>%
  select(site_id) %>%
  left_join(info, by = 'site_id') %>%
  select(site_id, ref_id,ref_pos, ref_allele, major_allele,
         minor_allele,
         locus_type, gene_id, snp_type, site_type, amino_acids) %>%
  # select(locus_type) %>% table
  HMVAR::determine_snp_effect() %>%
  HMVAR::determine_substitution_type()
info


