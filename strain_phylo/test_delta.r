#!/usr/bin/env Rscript
# setwd("/cashew/users/sur/exp/fraserv/2019/today2")
# All the delta code here is from:
# https://github.com/mrborges23/delta_statistic/blob/master/code.R

library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Calculate and test delta in a phylogenetic tree and trait vector."))
  
  # Positional arguments
  p <- add_argument(p, "tree_file",
                    help = paste("Newick file with phylogenetic tree."),
                    type = "character")
  
  # Optional arguments
  p <- add_argument(p, "--map_file",
                     help = paste("File with ID column and Group column containing",
                                  "trait. ID column must match tip labels of tree."),
                     type = "character",
                     default = "map.txt")
  p <- add_argument(p, "--nperm",
                    help = paste("Number of permutations to make."),
                    type = "numeric",
                    default = 100)
  p <- add_argument(p, "--output",
                    help = paste("Filename to write results. If <--type multi>, then",
                                 "this is the directory to write the output."),
                    type = "character",
                    default = "delta.txt")
  p <- add_argument(p, "--plot",
                    help = paste("Path to plot permutation results.",
                                 "If nothing passed no plot will be generated.",
                                 "If <--type multi> svg plots will always be",
                                 "generated."),
                    default = NULL,
                    type = "character")
  p <- add_argument(p, "--seed",
                    help = paste("Seed number for permutations"),
                    type = "numeric",
                    default = NULL)
  p <- add_argument(p, "--type",
                    help = paste("Either a single tree file ('single') or",
                                 "a directory with multiple tree files ('multi').",
                                 "Just one map file for all the trees is accepted."),
                    default = "single",
                    type = "character")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  if(!(args$type %in% c('single', 'multi'))){
    stop("ERROR: type must be 'single' or 'multi'.")
  }
  args$njobs <- 20
  
  return(args)
}

args <- process_arguments()
# args <- list(tree_file = "input_test/",
#              map_file = "map.txt",
#              nperm = 10,
#              output = "output/",
#              plot = "plot.svg",
#              seed = 1313,
#              type = "multi",
#              njobs = 20)

library(tidyverse)
library(ape)
library(expm)
library(rSubmitter)

################################################################################
#NENTROPY
#returns the node entropies by calculating sum of the state entropies
#prob: matrix of state probabilities

nentropy <- function(prob) {
  
  k              <- ncol(prob)                       #number of states
  prob[prob>1/k] <- prob[prob>1/k]/(1-k) - 1/(1-k)   #state entropies
  tent           <- apply(prob,1,sum)                #node entropy
  
  #correct absolute 0/1
  tent[tent == 0] <- tent[tent == 0] + runif(1,0,1)/10000
  tent[tent == 1] <- tent[tent == 1] - runif(1,0,1)/10000
  
  return(tent)
}

#FUNCTION FOR BAYESIAN INFERENCES
#bayesian inferences on the node entropies 
#l0: rate parameter of the exponential prior distribution
#se: standard deviation of the proposal distribution 
#a:  alpha parameter (beta likelihood)
#b:  beta paramter (beta likelihood)
#x:  node entropies

lpalpha <- function(a,b,x,l0) {          #log posterior alpha
  N  <- length(x)
  lp <- N*(lgamma(a+b)-lgamma(a)) - a*(l0-sum(log(x)))
  return(lp)
}

lpbeta  <- function(a,b,x,l0) {          #log posterior beta
  N  <- length(x)
  lp <- N*(lgamma(a+b)-lgamma(b)) - b*(l0-sum(log(1-x)))
  return(lp)
}

mhalpha <- function(a,b,x,l0,se) {       #metropolis hastings alpha
  a0 <- a
  a1 <- exp(rnorm(1,log(a0),se))
  
  r  <- min(1, exp(lpalpha(a1,b,x,l0) - lpalpha(a0,b,x,l0) ) )
  
  while (is.na(r) == T) {
    a1 <- exp(rnorm(1,log(a0),se))
    r  <- min(1, exp(lpalpha(a1,b,x,l0) - lpalpha(a0,b,x,l0) ) )
  }
  
  if (runif(1) < r) {
    return(a1) 
  } else {
    return(a0)
  }
}

mhbeta  <- function(a,b,x,l0,se) {      #metropolis hastings beta
  b0 <- b
  b1 <- exp(rnorm(1,log(b0),se))
  
  r  <- min(1, exp(lpbeta(a,b1,x,l0) - lpbeta(a,b0,x,l0) ) )
  
  while (is.na(r) == T) {
    b1 <- exp(rnorm(1,log(b0),se))
    r  <- min(1, exp(lpbeta(a,b1,x,l0) - lpbeta(a,b0,x,l0) ) )
  }  
  
  if (runif(1) < r) {
    return(b1)
  } else {
    return(b0)
  }
}

#MCMC
#Markov chain monte carlo scheme using the conditional posteriors of alpha and beta
#alpha: initial value of alpha
#beta: initial values of beta
#x: node entropies
#sim: number of iterations
#thin: controles the number of saved iterations = sim/thin
#burn: number of iterates to burn

emcmc <- function(alpha,beta,x,l0,se,sim,thin,burn) {
  
  usim <- seq(burn,sim,thin)
  gibbs <- matrix(NA,ncol=2,nrow=length(usim))
  p <- 1
  
  for (i in 1:sim) {
    alpha <- mhalpha(alpha,beta,x,l0,se)
    beta  <- mhbeta(alpha,beta,x,l0,se)
    
    if (i == usim[p]) {
      gibbs[p,] <- c(alpha,beta)
      p <- p+1
    }
  }  
  return(gibbs)
}

#RATE MATRIX FOR TRAIT EVOLUTION. K=2 TO 5
ratematrix <- function(pi,rho){
  
  k <- length(pi)
  
  if (k==2){
    r <- c(pi[1]*0     ,pi[2]*rho[1],
           pi[1]*rho[1],pi[2]*0)
  }
  
  if (k==3){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[3],
           pi[1]*rho[2],pi[2]*rho[3],pi[3]*0 )
  }
  
  if (k==4){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],pi[4]*rho[3],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[4],pi[4]*rho[5],
           pi[1]*rho[2],pi[2]*rho[4],pi[3]*0     ,pi[4]*rho[6],
           pi[1]*rho[3],pi[2]*rho[5],pi[3]*rho[6],pi[4]*0 )
  }  
  
  if (k==5){
    r <- c(pi[1]*0     ,pi[2]*rho[1],pi[3]*rho[2],pi[4]*rho[3] ,pi[5]*rho[4],
           pi[1]*rho[1],pi[2]*0     ,pi[3]*rho[5],pi[4]*rho[6] ,pi[5]*rho[7],
           pi[1]*rho[2],pi[2]*rho[5],pi[3]*0     ,pi[4]*rho[8] ,pi[5]*rho[9],
           pi[1]*rho[3],pi[2]*rho[6],pi[3]*rho[8],pi[4]*0      ,pi[5]*rho[10],
           pi[1]*rho[4],pi[2]*rho[7],pi[3]*rho[9],pi[4]*rho[10],pi[5]*0)
  }
  
  R <- matrix(r,ncol=k,nrow=k) 
  diag(R) <- -rowSums(R)
  
  return(R)
}

#RTRAIT
#simulates the evolution of a trait in a given tree
# tree: metric-tree
# R: rate matrix
# nstates: number of states

rtrait <- function(tree,R,nstates) {
  
  nspecis <- length(tree$tip.label)
  
  #tree
  edge <- cbind(tree$edge,tree$edge.length)
  
  ancestral <- rep(NA,2*nspecies-1) 
  ancestral[nspecies+1] <- sample(1:nstates,1,prob=pi) 
  
  #rate change
  inode <- nspecies+1
  while (sum(is.na(ancestral)) > 0) {
    
    inode1 <-  edge[which(edge[,1]==inode)[1],2]
    inode2 <-  edge[which(edge[,1]==inode)[2],2]
    bl1 <- edge[which(edge[,1]==inode)[1],3]
    bl2 <- edge[which(edge[,1]==inode)[2],3]
    
    astate <- rep(0,nstates)
    astate[ancestral[inode]] <- 1 
    
    ancestral[inode1] <- sample(1:nstates,1,prob=astate%*%expm(R*bl1))
    ancestral[inode2] <- sample(1:nstates,1,prob=astate%*%expm(R*bl2))
    
    inode <- inode+1
  }
  return(ancestral[1:nspecies])
  
}

#DELTA
#calculate delta statistic
#trait: trait vector 
delta <- function(trait, tree,lambda0,se,sim,thin,burn) {
  
  ar <- ace(trait,tree,type="discret",method="ML",model="ARD")$lik.anc
  x  <- nentropy(ar)
  mc1    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
  mc2    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
  mchain <- rbind(mc1,mc2)
  deltaA <- mean(mchain[,2]/mchain[,1])
  
  return(deltaA)
}

#' Title
#'
#' @param tre 
#' @param trait 
#' @param verbose 
#'
#' @return
#' @export
#' @author Sur
#'
#' @examples
calculate_delta <- function(tre, trait, verbose = FALSE){
  tre <- ape::multi2di(tre)
  tre$edge.length[ tre$edge.length <= 0 ] <- min(tre$edge.length[ tre$edge.length > 0])
  if(verbose)
    cat("\tCalculating delta...\n")
  tre.delta <- tryCatch(tre.delta <- delta(trait = trait, tree = tre,
                                           lambda0 = 0.1, se = 0.5,
                                           sim = 10000, thin = 10, burn = 200),
                        error = function(e){ NA })
  
}

#' Title
#'
#' @param tree_file 
#' @param trait 
#' @param seed 
#' @param nperm 
#' @param plot 
#' @param multi 
#' @param output 
#'
#' @return
#' @export
#' @author Sur
#'
#' @examples
delta_from_tree <- function(tree_file, trait, seed = NA, nperm = 100, plot = NA, multi = FALSE, output = "./"){
  # tree_file <- "a/b/c.d.tre"
  # tree_file <- "a/b/c"
  spec_name <- str_split(basename(tree_file), pattern = '[.]')[[1]]
  spec_name <- paste(spec_name[1:(length(spec_name) - 1)], collapse = ".")
  
  tre <- read.tree(tree_file)
  trait <- trait[ tre$tip.label ]
  if(any(is.na(trait))){
    stop("ERROR: not all tree tips in map")
  }
  
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  # Sample delta
  tre.delta <- calculate_delta(tre, trait = trait)
  
  # Permute
  delta.perms <- 1:nperm %>%
    purrr::map_dbl(~calculate_delta(tre = tre,
                                    trait = sample(trait, replace = FALSE),
                                    verbose = TRUE))
  
  # Plot
  if(multi){
    plot <- file.path(output, paste0(spec_name, ".svg"))
  }
  if(!is.na(plot)){
    p1 <- tibble(delta = delta.perms) %>%
      ggplot(aes(x = delta)) +
      geom_histogram(bins = 20) +
      geom_vline(xintercept = tre.delta,
                 col = "red", size = 2) +
      AMOR::theme_blackbox()
    ggsave(plot, p1, width = 6, height = 4)
  }
  
  
  
  
  # Get result
  Res <- tibble(name = spec_name,
                delta = tre.delta,
                pval = sum(delta.perms >= tre.delta) / sum(!is.na(delta.perms)),
                k = length(unique(trait)),
                N = length(trait),
                nperms = nperm,
                succ_perms = sum(!is.na(delta.perms)),
                seed = seed,
                mean_delta = mean(delta.perms),
                sd_delta = sd(delta.perms),
                median_delta = median(delta.perms),
                Q1 = quantile(delta.perms, 0.25),
                Q3 = quantile(delta.perms, 0.75),
                ci_bottom = quantile(delta.perms, 0.025),
                ci_top = quantile(delta.perms, 0.975))
  
  return(Res)
}
################################################################################

# Create trait vector
trait <- read_tsv(args$map_file)
trait <- setNames(trait$Group, trait$ID)

if(args$type == "single"){
  # args <- list(tree_file = "Bacteroides_vulgatus_57955.tre",
  #              map_file = "map.txt",
  #              nperm = 10,
  #              output = "output/",
  #              plot = "plot.svg",
  #              seed = 1313,
  #              type = "single",
  #              njobs = 20)
  
  Res <- delta_from_tree(tree_file = args$tree_file,
                         trait = trait, seed = args$seed,
                         nperm = args$nperm, plot = args$plot,
                         multi = FALSE)
  write_tsv(Res, args$output)
}else if(args$type == "multi"){
  # args <- list(tree_file = "input_test/",
  #              map_file = "map.txt",
  #              nperm = 10,
  #              output = "output/",
  #              plot = "plot.svg",
  #              seed = 1313,
  #              type = "multi",
  #              njobs = 20)
  
  files <- list.files(args$tree_file, full.names = TRUE, recursive = FALSE)
  if(length(files) < 1){
    stop("ERROR: no files in input directory", call. = TRUE)
  }
  
  if(!dir.exists(args$output)){
    dir.create(args$output)
  }
  
  Res <- superApply(x = files, FUN = delta_from_tree,
                    trait = trait, seed = args$seed,
                    nperm = args$nperm, multi = TRUE,
                    output = args$output,
                    tasks = args$njobs,
                    workingDir = "workdir/",
                    clean = FALSE,
                    extraBashLines = "module load R/3.6.1")
  filename <- file.path(args$output, "delta_test.txt")
  write_tsv(Res %>% bind_rows(), filename)
}else{
  stop("ERROR: bad --type", call. = TRUE)
}