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


setwd("/home/sur/micropopgen/exp/2019/today")
library(tidyverse)
library(HMVAR)



args <- list(depth_thres = 1,
             freq_thres =  0.5,
             min_cov =  0.8,
             keep_last_codon = TRUE,
             outdir = 'output/',
             genomes_dir = "/home/sur/micropopgen/data/genomes/midas_db_v1.2/hmp.subsite/",
             map_file = "midas/map.txt",
             midas_dir = "midas/merged.snps/Veillonella_parvula_57794/",
             missing_as = 'gap')

if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}
map <- read_tsv(args$map_file,
                col_types = cols(.default = col_character())) %>%
  select(sample = ID, Group)


map

for(midas_dir in args$midas_dir){
  midas_dir <- args$midas_dir[1]
  spec <- basename(midas_dir)
  cat(spec, "\n")
  
  # Read data
  Dat <- read_data(spec = spec, midas_dir = midas_dir, genomes_dir = args$genomes_dir)
  # map %>% filter(sample %in% colnames(Dat$midas$freq)) %>% select(Group) %>% table
  
  if(ncol(Dat$midas$freq) > 2){
    # At least two samples for analysis
    
    # For each gene, extract MSA
    Res <- Dat$midas$info %>%
      # head(1000) %>%
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
        f <- Dat$freq %>% dplyr::filter(site_id %in% i$site_id)
        d <- Dat$depth %>% dplyr::filter(site_id %in% i$site_id)
        
        # Calculate coverage and discard samples with uncovered gene
        gene_cov <- colSums((d %>% select(-site_id)) > depth_thres) / nrow(i)
        if(sum(gene_cov >= min_cov) < min_cov*length(gene_cov)){
          gene_cov <- tibble(gene = gene, sample = names(gene_cov), coverage = gene_cov)
          return(list(Coverage = gene_cov, aln = NULL))
        }
        d <- d %>% select(site_id, names(which(gene_cov >= min_cov)))
        f <- f %>% select(site_id, names(which(gene_cov >= min_cov)))
        gene_cov <- tibble(gene = gene, sample = names(gene_cov), coverage = gene_cov)
        
        # Match
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
        
        # p1 <- ggplot(dat, aes(x = Group, y = freq)) +
        #   facet_wrap(~site_id, scales = "free_y") +
        #   geom_point(position = position_jitter(width = 0.3))
        # p1
        
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
        # outfile <- paste0(gene, ".aln")
        # outfile <- file.path(outdir, outfile)
        # seqinr::translate(aln$seq$SRS143876)
        # seqinr::write.fasta(aln$seq, aln$nam, file.out = outfile)
        
        return(list(Coverage = gene_cov, aln=aln$seq))
      }, freq=Dat$midas$freq, depth=Dat$midas$depth, map=map,
      depth_thres=args$depth_thres, freq_thres=args$freq_thres, min_cov=args$min_cov,
      genome_feats=Dat$genome_feats, missing_as=args$missing_as,
      genome_fasta=Dat$genome_fasta)
    
  }else{
    Res <- NULL
  }
  
  
}


# 
# cov <- Res %>%
#   map_dfr(~ .$Coverage) %>%
#   spread(sample, coverage)
# genes <- cov$gene
# cov <- cov %>%
#   select(-gene) %>%
#   as.matrix
# row.names(cov) <- genes
# cov <- cov[ rowSums(cov >= min_cov) >= (ncol(cov) * min_cov), ]
# cov <- cov[ , colSums(cov >= min_cov) >= (nrow(cov) * min_cov) ]
# 
# map %>%
#   filter(sample %in% colnames(cov)) %>%
#   select(Group) %>%
#   table
# 
# 
# Alns <- Res[row.names(cov)] %>%
#   map(function(l, samples){
#     aln <- l$aln
#     aln <- aln[intersect(names(aln), samples)]
#     return(aln)
#   }, samples = colnames(cov))
# 
# dir.create(outdir)
# names(Alns) %>%
#   map(function(gene, Alns, outdir="./"){
#     # cat(gene, "\n")
#     outfile <- paste0(gene, ".aln")
#     outfile <- file.path(outdir, outfile)
#     