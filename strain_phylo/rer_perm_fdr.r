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

library(tidyverse)

# args <- list(permdir = "perm_res/rertest/Bacteroides_vulgatus_57955/",
#              res = "../2019-11-19.rertest/gut/rer/rer_output/rertest/Bacteroides_vulgatus_57955/Bacteroides_vulgatus_57955.cors.txt",
#              output = "Bacterioides_vulgatus_57955.rer.fdr.txt")
opts <- commandArgs(trailingOnly=TRUE)
args <- list(permdir = opts[1],
             res = opts[2],
             output = opts[3])

# Read data
cat("Read original results...\n")
Res <- read_tsv(args$res,
                col_types = cols(gene_id = col_character()))
cat("Read permutations...\n")
Perms <- list.files(args$permdir, full.names = TRUE) %>%
  map(~read_tsv(., col_types = cols(gene_id = col_character()))$P) %>%
  bind_cols %>%
  bind_cols(orig = Res$P)

# Calculate FDR
cat("Calculating FDR...\n")
Res$FDR <- NA
for(i in 1:nrow(Res)){
  pval <- Res$P[i]

  fdr <- Perms %>%
    map_dbl(function(x){
      x <- x[!is.na(x)]
      sum(x <= pval) / length(x)
    }) %>% mean
  Res$FDR[i] <- fdr
}

# Write output
cat("Writing output...\n")
write_tsv(Res, path = args$output)
