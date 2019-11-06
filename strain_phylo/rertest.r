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

##############################################
library(argparser)

process_arguments <- function(){
  p <- arg_parser(paste("Code that runs RERconverge test on a binary phenotype"))
  
  # Positional arguments
  p <- add_argument(p, "indir",
                    help = paste("Directory containing one all trees from a given species",
                                 "or a list of subdirectories (one per species) with each",
                                 "one containing all tree species of the corresponding",
                                 "species. CURRENTLY ONLY FIRST OPTION AVAILABLE."),
                    type = "character")
  p <- add_argument(p, "master_tree",
                    help = paste("Path to master phylogenetic tree file. Must include",
                                 "all species"),
                    type = "character")
  
  # Optional arguments
  # p <- add_argument(p, "--type",
  #                    help = paste("Indicates if <indir> corresponds to a single species",
  #                                 "(single), or to a directory with subdirectories for",
  #                                 "multiples species (multi)."),
  #                    type = "character",
  #                    default = "single")
  p <- add_argument(p, "--map_file",
                    help = paste("Path to file with map. Should have 'ID' and 'group'",
                                 "columns."),
                    type = "character",
                    default = "map.txt")
  p <- add_argument(p, "--outdir",
                    help = paste("Directory path to store outputs."),
                    default = "output/",
                    type = "character")
  p <- add_argument(p, "--focal_phenotype",
                    help = paste("Phenotypic value to be tested. Internally, tree tips",
                                 "with this value will be set to one, and to zero",
                                 "otherwise. If nothing passed, the script assumes that",
                                 "there are only two values in the Group column of the",
                                 "<--map_file>."),
                    default = NULL,
                    type = "character")
                     
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  # if(!(args$type %in% c('single', 'multi'))){
  #   stop("ERROR: --type must be 'single' or 'multi'", call. = TRUE)
  # }
  
  # Set other arguments
  args$scale <- TRUE
  args$weight <- TRUE
  args$type <- "single"
  args$transform <- 'sqrt'
  args$clade <- "all"
  args$min.sp <- 5
  args$min.pos <- 2
  
  return(args)
}

#' Computes the association statistics between RER from
#' \code{\link{getAllResiduals}} and a phenotype paths vector
#' made with \code{\link{tree2Paths}} or \code{\link{char2Paths}}.
#' Modified by Sur to handle errors.
#' @param RERmat RER matrix returned by \code{\link{getAllResiduals}}
#' @param charP phenotype vector returned by \code{\link{tree2Paths}}
#' or \code{\link{char2Paths}}
#' @param method Method used to compute correlations. Accepts the
#' same arguments as \code{\link{cor}}. Set to "auto" to select automatically
#' based on the number of unique values in charP. This will also auto set
#' the winsorization for Pearson correlation. Set winsorizetrait=some number
#' and winsorizeRER=some number to override
#' @param min.sp Minimum number of species that must be present for a gene
#' @param min.pos Minimum number of species that must be present in the
#' foreground (non-zero phenotype values)
#' @param winsorizeRER Winsorize RER values before computing Pearson correlation.
#' winsorizeRER=3 will set the 3 most extreme RER values at each end of each row to the value closest to 0.
#' @param winsorizetrait Winsorize trait values before computing Pearson
#' correlation. winsorizetrait=3 will set the 3 most extreme trait values at each
#' end to the value closest to 0.
#' @param weighted perform weighted correlation. This option needs to be set if
#' the clade weights computed in \code{\link{foreground2Tree}(wholeClade=T)}
#' are to be used. This setting will treat the clade a single observation for
#' the purpose of p-value estimation.
#' @note  winsorize is in terms of number of observations at each end, NOT quantiles
#' @return A list object with correlation values, p-values, and the number of
#' data points used for each tree
#' @export
getAllCor <- function(RERmat, charP, method="auto",min.sp=10, min.pos=2,
                      winsorizeRER=NULL, winsorizetrait=NULL, weighted=F){
  # RERmat <- rerw
  # charP <- phenv
  # method <- 'k'
  # min.sp <- 5
  # min.pos <- 2
  # winsorizeRER <- NULL
  # winsorizetrait <- NULL
  # weighted <- TRUE
  
  RERna=(apply(is.na(RERmat),2,all))
  iicharPna=which(is.na(charP))
  if(!all(RERna[iicharPna])){
    warning("Species in phenotype vector are a subset of the those used for RER computation. For best results run getAllResiduals with the useSpecies")
  }
  if (method=="auto"){
    lu=length(unique(charP))
    if(lu==2){
      method="k"
      message("Setting method to Kendall")
    }
    else if (lu<=5){
      method="s"
      message("Setting method to Spearman")
    }
    else{
      method="p"
      message("Setting method to Pearson")
      if(is.null(winsorizeRER)){
        message("Setting winsorizeRER=3")
        winsorizeRER=3
      }
      if(is.null(winsorizetrait)){
        message("Setting winsorizetrait=3")
        winsorizetrait=3
      }
    }
  }
  win=function(x,w){
    xs=sort(x[!is.na(x)], decreasing = T)
    xmax=xs[w]
    xmin=xs[length(xs)-w+1]
    
    x[x>xmax]=xmax
    x[x<xmin]=xmin
    x
  }
  corout=matrix(nrow=nrow(RERmat), ncol=3)
  rownames(corout)=rownames(RERmat)
  
  colnames(corout)=c("Rho", "N", "P")
  
  for( i in 1:nrow(corout)){
    cat(row.names(corout)[i], "\n")
    
    if(((nb<-sum(ii<-(!is.na(charP)&!is.na(RERmat[i,]))))>=min.sp)){
      if (method!="p"&&sum(charP[ii]!=0)<min.pos){
        next
      }
      
      if(!weighted){
        
        x=RERmat[i,]
        
        #winsorize
        indstouse=which(!is.na(x) & !is.na(charP))
        if(!is.null(winsorizeRER)){
          x=win(x[indstouse], winsorizeRER)
        }else{
          x=x[indstouse]
        }
        if(!is.null(winsorizetrait)){
          y=win(charP[indstouse], winsorizetrait)
        }else{
          y=charP[indstouse]
        }
        
        cres=cor.test(x, y, method=method, exact=F)
        corout[i,1:3]=c(cres$estimate, nb, cres$p.value)
      }
      else{
        charPb=(charP[ii]>0)+1-1
        
        weights=charP[ii]
        weights[weights==0]=1
        
        cres <- tryCatch(wtd.cor(RERmat[i,ii], charPb, weight = weights, mean1 = F),
                         error = function(e){
                           cat("\tskipping...\n")
                         })
        if(is.null(cres)){
          next
        }else{
          corout[i, 1:3]=c(cres[1], nb, cres[4])
        }
      }
    }
    else{
      #show(i)
      #show(c(nb, charP[ii]))
    }
    
  }
  
  corout=as.data.frame(corout)
  corout$p.adj=p.adjust(corout$P, method="BH")
  corout
}
##############################################

args <- process_arguments()

library(tidyverse)
library(RERconverge)


# tree_dir <- "output/gene_trees/"
# tree_tab_file <- "trees_tab.txt"
# map_file <- "../../../../data/gathered_results/2019a.gut/map.txt"

# Read simple data
master_tre <- ape::read.tree(args$master_tree)

map <- read_tsv(args$map_file)
map <- setNames(map$Group, map$ID)
map <- map[ Trees$masterTree$tip.label ]
if(is.na(args$focal_phenotype)){
  args$focal_phenotype <- map[1]""
}
map <- 1*(map == args$foca_phenotype)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(outdir)
}

for(specdir in args$indir){
  spec <- basename(specdir)
  cat(spec, "\n")
  
  # Read trees
  Trees <- readTrees(tree_tab_file, masterTree = master_tre)
  filename <- file.path(args$oudir, "Trees.dat")
  save(Trees, file = filename)
  
  op <- par()
  # Calculate RERs
  rerw <- getAllResiduals(Trees, transform = args$transform, weighted = args$weight, scale = args$scale, plot = FALSE)
  filename <- file.path(args$oudir, "rerw.dat")
  save(rerw, file = filename)
  
  # Prepare binary tree
  bintre <- foreground2Tree(foreground = names(map)[map == 1],
                            treesObj = Trees,
                            plotTree = FALSE,
                            clade = args$clade,
                            weighted = TRUE)
  
  # Prepare phenotypes
  phenv <- tree2Paths(tree = bintre, treesObj = Trees)
  
  # Calculate correlatiob between phenotypes and rers
  cor.res <- correlateWithBinaryPhenotype(RERmat = rerw, charP = phenv, min.sp = args$min.sp, min.pos = args$min.pos)
  # cor.res <- getAllCor(rerw, phenv, 5, 2, method = "k", weighted=TRUE)
  cor.res <- cor.res[ order(cor.res$P), ]
  filename <- file.path(args$oudir, "cors.txt")
  write_tsv(cor.res, path = filename)
}




