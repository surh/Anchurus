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
  p <- add_argument(p, "trees",
                    help = paste("Table file containing all trees from a given species",
                                 "or a directory containing multiple files (one per",
                                 "species) with eachcone containing all tree species",
                                 "of the corresponding species. CURRENTLY ONLY FIRST",
                                 "OPTION AVAILABLE."),
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
                    help = paste("Path to file with map. Should have 'ID' and 'Group'",
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
  p <- add_argument(p, "--spec",
                    help = paste("Species name. If passed it will be used to construct",
                                 "filenames."),
                    type = "character",
                    default = NULL)
  
  # Read arguments
  cat("Processing arguments...\n")
  args <- parse_args(p)
  
  # Process arguments
  # if(!(args$type %in% c('single', 'multi'))){
  #   stop("ERROR: --type must be 'single' or 'multi'", call. = TRUE)
  # }
  if(is.na(args$spec)){
    args$spec <- NULL
  }
  
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

getAllResiduals <- function(treesObj, cutoff=NULL, transform="sqrt", weighted=T,
                            useSpecies=NULL,  min.sp=10, scale=T, 
                            doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05, plot=T){
  
  # treesObj <- Trees
  # cutoff <- NULL
  # transform <- args$transform
  # weighted <- args$weight
  # useSpecies <- NULL
  # min.sp <- 10 
  # scale <- args$scale
  # doOnly <- NULL
  # maxT <- NULL
  # scaleForPproj <- FALSE
  # mean.trim <- 0.05
  # plot <- FALSE
  # 
  
  if(is.null(cutoff)){
    cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
    message(paste("cutoff is set to", cutoff))
  }
  if (weighted){
    weights=RERconverge:::computeWeightsAllVar(treesObj$paths, transform=transform, plot=plot)
    residfunc=RERconverge:::fastLmResidMatWeighted
  }else{
    residfunc=fastLmResidMat
  }
  # residfunc=naresid
  
  if (is.null(useSpecies)){
    useSpecies=treesObj$masterTree$tip.label
    #mappedEdges=trees$mappedEdges
  }
  if(is.null(maxT)){
    maxT=treesObj$numTrees
  }
  if(transform!="none"){
    transform=match.arg(transform,c("sqrt", "log"))
    transform=get(transform)
  }else{
    transform=NULL
  }
  
  
  
  #cm is the names of species that are included in useSpecies and the master tree
  cm=intersect(treesObj$masterTree$tip.label, useSpecies)
  sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
  if (length(sp.miss) > 0) {
    message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
                                                                                 collapse = ",")))
  }
  
  rr=matrix(nrow=nrow(treesObj$paths), ncol=ncol(treesObj$paths))
  
  #maximum number of present species
  maxn=rowSums(treesObj$report[,cm])
  
  if(is.null(doOnly)){
    doOnly=1
  }else{
    maxT=1
  }
  skipped=double(nrow(rr))
  skipped[]=0
  
  for (i in doOnly:(doOnly+maxT-1)){
    
    if(sum(!is.na(rr[i,]))==0&&!skipped[i]==1){
      
      
      #get the ith tree
      tree1=treesObj$trees[[i]]
      
      #get the common species, prune and unroot
      both=intersect(tree1$tip.label, cm)
      if(length(both)<min.sp){
        next
      }
      tree1=unroot(pruneTree(tree1,both))
      
      #do the same for the refTree
      
      
      #find all the genes that contain all of the species in tree1
      allreport=treesObj$report[,both]
      ss=rowSums(allreport)
      iiboth=which(ss==length(both)) #this needs to be >1
      if (length(iiboth) < 2) {
        message(paste("Skipping i =",i,"(no other genes with same species set)"))
        next
      }
      
      nb=length(both)
      ai=which(maxn[iiboth]==nb)
      
      
      message(paste("i=", i))
      
      
      if(T){
        
        ee=RERconverge:::edgeIndexRelativeMaster(tree1, treesObj$masterTree)
        
        ii= treesObj$matIndex[ee[, c(2,1)]]
        
        allbranch=treesObj$paths[iiboth,ii]
        if (is.null(dim(allbranch))) {
          message(paste("Issue with gettiing paths for genes with same species as tree",i))
          return(list("iiboth"=iiboth,"ii"=ii))
        }
        
        if(weighted){
          allbranchw=weights[iiboth,ii]
        }
        if(scaleForPproj){
          nv=apply(scaleMatMean(allbranch), 2, mean, na.rm=T, trim=mean.trim)
        }else{
          nv=apply(allbranch, 2, mean, na.rm=T, trim=mean.trim)
        }
        
        iibad=which(allbranch<cutoff)
        #don't scale
        #allbranch=scaleMat(allbranch)
        if(!is.null(transform)){
          nv=transform(nv)
          allbranch=transform(allbranch)
        }
        allbranch[iibad]=NA
        
        
        
        # cat("\thello...\n")
        if(!scale){
          if(!weighted){
            proj=residfunc(allbranch[ai, ,drop=F], model.matrix(~1+nv))
            
          }else{
            
            proj=residfunc(allbranch[ai, ,drop=F], model.matrix(~1+nv), allbranchw[ai, ,drop=F])
            
          }
        }else{
          if(!weighted){
            proj=residfunc(allbranch[, ,drop=F], model.matrix(~1+nv))
          }else{
            # proj=residfunc(allbranch[, ,drop=F], model.matrix(~1+nv),allbranchw)
            
            proj <- tryCatch(residfunc(allbranch[, ,drop=F], model.matrix(~1+nv),allbranchw),
                             error = function(e){
                               cat("\tskipping...\n")
                             })
            if(is.null(proj)){
              next
            }
          }
          proj=scale(proj, center = F)[ai, , drop=F]
          
        }
        
        
        #we have the projection
        
        
        
        rr[iiboth[ai],ii]=proj
        
      }
      
    }}
  message("Naming rows and columns of RER matrix")
  rownames(rr)=names(treesObj$trees)
  colnames(rr)=namePathsWSpecies(treesObj$masterTree)
  rr
}


#' @param file The path to the tree file
#' @param  max.read This function takes a while for a whole genome, so max.read is useful for testing
#' @param  masterTree (optional) User can specify a master tree; only the topology will be used, and branch lengths will be inferred from gene trees.
#' @return A trees object of class "treeObj"
#' @param  masterTree (optional) User can specify a master tree. Recommended only when
#' the number of available gene trees with all species is small.
#' @param  minTreesAll The minimum number of trees with all species present in order to estimate
#' master tree edge lengths (default 20).
#' @return A trees object of class "treeObj"
#' @export
#' @author RERconverge package. Modified by Sur
readTrees=function(file, max.read=NA, masterTree=NULL, minTreesAll=20){
  
  tmp=scan(file, sep="\t", what="character", quiet = T)
  message(paste0("Read ",length(tmp)/2, " items", collapse=""))
  trees=vector(mode = "list", length = min(length(tmp)/2,max.read, na.rm = T))
  treenames=character()
  maxsp=0; # maximum number of species
  allnames=NA # unique tip labels in gene trees
  
  #create trees object, get species names and max number of species
  for ( i in 1:min(length(tmp),max.read*2, na.rm = T)){
    if (i %% 2==1){
      treenames=c(treenames, tmp[i])
    }
    else{
      trees[[i/2]]=tryCatch(unroot(read.tree(text=tmp[i])),
                            error = function(e) {
                              message('Cannot parse tree for the following gene: ',treenames[i/2]);
                              stop()
                            })
      #reduce to species present in master tree
      if (!is.null(masterTree)) {
        trees[[i/2]] = pruneTree(trees[[i/2]],intersect(trees[[i/2]]$tip.label,masterTree$tip.label))
      }
      
      #check if it has more species
      # cat("===========\n")
      # cat("maxsp:", maxsp, "\n")
      # cat("len(allnames):", length(allnames), "\n")
      if(length(trees[[i/2]]$tip.label)>maxsp){
        maxsp=length(trees[[i/2]]$tip.label)
        allnames=trees[[i/2]]$tip.label
      }
      #check if it has new species
      if (sum(trees[[i/2]]$tip.label %in% allnames == F) > 0) {
        allnames = unique(c(allnames,trees[[i/2]]$tip.label))
        maxsp = length(allnames) - 1
        
      }
      #if(length(trees[[i/2]]$tip.label)>maxsp){
      #  maxsp=length(trees[[i/2]]$tip.label)
      #  allnames=trees[[i/2]]$tip.label
      #}
    }
    
  }
  allnames = allnames[!is.na(allnames)]
  
  maxsp <- length(allnames)
  
  names(trees)=treenames
  treesObj=vector(mode = "list")
  treesObj$trees=trees
  treesObj$numTrees=length(trees)
  treesObj$maxSp=maxsp
  
  message(paste("max is", maxsp))
  
  report=matrix(nrow=treesObj$numTrees, ncol=maxsp)
  colnames(report)=allnames
  
  rownames(report)=treenames
  for ( i in 1:nrow(report)){
    ii=match(allnames, trees[[i]]$tip.label)
    report[i,]=1-is.na(ii)
    
  }
  treesObj$report=report
  
  
  
  ii=which(rowSums(report)==maxsp)
  
  #Create a master tree with no edge lengths
  if (is.null(masterTree)) {
    master=trees[[ii[1]]]
    master$edge.length[]=1
    treesObj$masterTree=master
  } else {
    
    master=pruneTree(masterTree, intersect(masterTree$tip.label,allnames))
    #prune tree to just the species names in the largest gene tree
    master$edge.length[]=1
    
    master=unroot(pruneTree(masterTree, intersect(masterTree$tip.label,allnames)))
    #prune tree to just the species names in the gene trees
    #master$edge.length[]=1
    
    treesObj$masterTree=master
  }
  
  
  
  
  treesObj$masterTree=rotateConstr(treesObj$masterTree, sort(treesObj$masterTree$tip.label))
  #this gets the abolute alphabetically constrained order when all branches
  #are present
  tiporder=RERconverge:::treeTraverse(treesObj$masterTree)
  
  #treesObj$masterTree=CanonicalForm(treesObj$masterTree)
  message("Rotating trees")
  
  for ( i in 1:treesObj$numTrees){
    
    treesObj$trees[[i]]=rotateConstr(treesObj$trees[[i]], tiporder)
    
  }
  
  
  
  ap=RERconverge:::allPaths(master)
  treesObj$ap=ap
  matAnc=(ap$matIndex>0)+1-1
  matAnc[is.na(matAnc)]=0
  
  paths=matrix(nrow=treesObj$numTrees, ncol=length(ap$dist))
  for( i in 1:treesObj$numTrees){
    #Make paths all NA if tree topology is discordant
    paths[i,]=tryCatch(RERconverge:::allPathMasterRelative(treesObj$trees[[i]], master, ap), error=function(err) NA)
    #calls matchAllNodes -> matchNodesInject
  }
  paths=paths+min(paths[paths>0], na.rm=T)
  treesObj$paths=paths
  treesObj$matAnc=matAnc
  treesObj$matIndex=ap$matIndex
  treesObj$lengths=unlist(lapply(treesObj$trees, function(x){sqrt(sum(x$edge.length^2))}))
  
  #require all species and tree compatibility
  #ii=which(rowSums(report)==maxsp)
  ii=intersect(which(rowSums(report)==maxsp),which(is.na(paths[,1])==FALSE))
  
  if (is.null(masterTree)) {
    if(length(ii)>=minTreesAll){
      message (paste0("estimating master tree branch lengths from ", length(ii), " genes"))
      tmp=lapply( treesObj$trees[ii], function(x){x$edge.length})
      
      allEdge=matrix(unlist(tmp), ncol=2*maxsp-3, byrow = T)
      allEdge=RERconverge:::scaleMat(allEdge)
      allEdgeM=apply(allEdge,2,mean)
      treesObj$masterTree$edge.length=allEdgeM
    }else {
      message("Not enough genes with all species present: master tree has no edge.lengths")
    }
  } else {
    message("Using user-specified master tree")
    
  }
  
  message("Naming columns of paths matrix")
  colnames(treesObj$paths)=namePathsWSpecies(treesObj$masterTree)
  class(treesObj)=append(class(treesObj), "treesObj")
  treesObj
}
##############################################

args <- process_arguments()
# args <- list(trees = "trees_tab.txt",
#              master_tree = "master_tree.tre",
#              map_file = "map.txt",
#              outdir = "output/",
#              focal_phenotype = "USA",
#              spec = "Bacteroides_ovatus_58035")
# # Set other arguments
# args$scale <- TRUE
# args$weight <- TRUE
# args$type <- "single"
# args$transform <- 'sqrt'
# args$clade <- "all"
# args$min.sp <- 5
# args$min.pos <- 2

library(tidyverse)
library(RERconverge)

# Read simple data
master_tre <- ape::read.tree(args$master_tree)

map <- read_tsv(args$map_file)
map <- setNames(map$Group, map$ID)

# Prepare output dir
if(!dir.exists(args$outdir)){
  dir.create(args$outdir)
}

for(specfile in args$trees){
  specfile <- args$trees[1]
  cat(basename(specfile), "\n")
  
  # Read trees
  Trees <- readTrees(specfile, masterTree = master_tre)
  filename <- file.path(args$outdir, paste(c(args$spec, "Trees.dat"), collapse = "."))
  cat("\twriting ", filename, "\n")
  save(Trees, file = filename)
  # load(file = filename)
  
  # Process map (probably need to change variable name for multi dir).
  map <- map[ Trees$masterTree$tip.label ]
  if(is.na(args$focal_phenotype)){
    args$focal_phenotype <- map[1]
  }
  map <- 1*(map == args$focal_phenotype)
  
  op <- par()
  # Calculate RERs
  rerw <- getAllResiduals(Trees, transform = args$transform, weighted = args$weight, scale = args$scale, plot = FALSE)
  filename <- file.path(args$outdir, paste(c(args$spec, "rerw.dat"), collapse = "."))
  cat("\twriting ", filename, "\n")
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
  # cor.res <- correlateWithBinaryPhenotype(RERmat = rerw, charP = phenv, min.sp = args$min.sp, min.pos = args$min.pos)
  cor.res <- getAllCor(RERmat = rerw, charP = phenv, min.sp = args$min.sp,
                       min.pos = args$min.pos, method = "k", weighted=TRUE)
  cor.res <- cor.res[ order(cor.res$P), ] %>% as_tibble(rownames = "gene_id")
  filename <- file.path(args$outdir, paste(c(args$spec, "cors.txt"), collapse = "."))
  cat("\twriting ", filename, "\n")
  write_tsv(cor.res, path = filename)
}




