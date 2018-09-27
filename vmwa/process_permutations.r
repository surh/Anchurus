#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(argparser)
###################################################
#' Process command line arguments
#'
#' @return List of arguments
#' 
#' @author Sur Herrera Paredes
#' 
#' @importFrom argparser arg_parser add_argument parse_args
#' @export
process_arguments <- function(){
  p <- argparser::arg_parser(paste0("Compares association results with permutations"))
  
  # Positional arguments
  p <- argparser::add_argument(p, "original", help = paste0("Results from real phenotypes."),
                               type = "character")
  p <- argparser::add_argument(p, "--perms", help = paste0("Files from permuted phenotypes."),
                               type = "character", nargs = Inf)
  
  # Read arguments
  args <- argparser::parse_args(p)
  
  return(args)
}
###################
# [d9/fe9ce6] Submitted process > permute_samples
# [e3/2105fb] Submitted process > run_vmwa
# [44/5e010f] Submitted process > run_vmwa_perms (3)
# [ab/a123bd] Submitted process > run_vmwa_perms (4)
# [6c/37c0df] Submitted process > run_vmwa_perms (5)
# [91/0fbd7d] Submitted process > run_vmwa_perms (2)
# [a0/2a0fdc] Submitted process > run_vmwa_perms (1)
# [cc/cef48a] Submitted process > run_vmwa_perms (8)
# [e8/c90fb5] Submitted process > run_vmwa_perms (9)
# [5d/e59161] Submitted process > run_vmwa_perms (7)
# [cf/ac6198] Submitted process > run_vmwa_perms (6)
# [d4/50c3fb] Submitted process > run_vmwa_perms (10)

# args <- list(original = "work/e3/2105fbf27d6251e3f1fa10eaaca0e3/results.txt",
#              perms = c("work/a0/2a0fdc992e27a091c7746640205318/results.txt",
#                        "work/91/0fbd7df0acb724a5d7d8de7210554d/results.txt",
#                        "work/44/5e010f29e6f37201c557ef29f9ef9f/results.txt",
#                        "work/ab/a123bdaa47b98b98b1875590946f4b/results.txt",
#                        "work/6c/37c0dfc52d52a8c28b24b2ac5ed63c/results.txt",
#                        "work/cf/ac6198716d364021810ceadd0fa7c5/results.txt",
#                        "work/5d/e5916139ea70e7949c7f59068fd848/results.txt",
#                        "work/cc/cef48adfb34eeb33c1a0692867fc98/results.txt",
#                        "work/e8/c90fb57f873140e325e3677ce0de26/results.txt",
#                        "work/d4/50c3fb4fdefe12ac738476e3f5917e/results.txt"))

args <- process_arguments()
# Read original
Res <- read.table(args$original, header = TRUE, stringsAsFactors = FALSE)
Res <- Res %>% mutate(count = 1)

# Compare to permutations
for(file in args$perms){
  Perm <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  row.names(Perm) <- Perm$SNP
  Perm <- Perm[ Res$SNP, ]
  Res$count <- Res$count + (abs(Perm$beta) >= abs(Res$beta))
}

# Calculate permutation P
Res <- Res %>% mutate(P = count / (length(args$perms) + 1))

# Sort
Res <- Res[ order(Res$P, decreasing = FALSE), ]

# Write output
p1 <- ggplot(Res, aes(x = P)) +
  geom_histogram(bins = 20) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("p.value_histogram.svg", p1, width = 6, height = 4)
write.table(Res, "results.txt", sep = "\t", quote = FALSE,
            col.names = TRUE, row.names = TRUE)
