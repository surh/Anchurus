library(tidyverse)
find_genes_and_samples <- function(cov, min_cov = 0.8){
  cat("\tSelecting genes...\n")
  genes <- cov$gene
  cov <- cov %>%
    select(-gene) %>%
    as.matrix
  row.names(cov) <- genes
  cov <- cov[ rowSums(cov >= min_cov) >= (ncol(cov) * min_cov), , drop = FALSE ]
  cov <- cov[ , colSums(cov >= min_cov) >= (nrow(cov) * min_cov), drop = FALSE ]
  
  list(samples = colnames(cov),
       genes = row.names(cov))
}

args <- list(min_cov = 0.8)

cov <- read_tsv("gene_covergages_gut/Bacteroides_vulgatus_57955.gene_coverage.txt")
# cov <- cov[1:50, 1:50]

dat <- cov %>%
  pivot_longer(-gene, names_to = "sample", values_to = "coverage")
dat

selected <- find_genes_and_samples(cov, min_cov = args$min_cov)

p1 <- dat %>%
  mutate(selected_sample = c("excluded", "selected")[1*(sample %in% selected$samples) + 1],
         selected_gene = c("excluded", "selected")[1*(gene %in% selected$genes) + 1]) %>%
  ggplot(aes(x = gene, y = sample)) +
  geom_tile(aes(fill = coverage)) +
  scale_fill_gradient(low = "white", high = "darkred") +
  facet_grid(selected_sample ~ selected_gene, scales = "free", space = "free") +
  theme(axis.text = element_blank())
# p1
ggsave("heatmap_Bvulgatus.png", p1, width = 10, height = 5)


# cov %>%
#   select(-gene) %>%
#   as.matrix() %>%
#   t() %>%
#   dist("manhattan") %>%
#   as.vector() %>%
#   summary()


























