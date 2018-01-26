library(AMOR)

# covfile <- commandArgs(trailingOnly = TRUE)[1]
# mapfile <- commandArgs(trailingOnly = TRUE)[2]
# covfile <- "coverage.txt"
# mapfile <- "map_for_comp.txt"
covfile <- opts[1]
mapfile <- opts[2]

# tab <- read.am("~/micropopgen/exp/2017/2017-08-07.abundances/merge.species/coverage.txt")
tab <- read.am(covfile)
tab

# map <- read.table("full_map.txt", header = TRUE, sep ="\t")
map <- read.table(mapfile, header = TRUE, sep ="\t")
head(map)
row.names(map) <- as.character(map$ID)

Dat <- create_dataset(Tab = tab$Tab, Map = droplevels(map[ samples(tab), ]))
Dat

comparisons <- combn(x = levels(Dat$Map$Group),m = 2)
comparisons

Res <- apply(comparisons, 2, function(x){
  #x <- comparisons[,1]
  
  cat(x, "\n")
  
  dat <- subset(Dat,Group %in% x)
  dat <- clean(dat)
  dat <- measurable_taxa(dat,min_reads_otu = 3, min_samples_otu = 4)
  thres <- 3
  
  Res <- NULL
  for(s in taxa(dat)){
    #s <- taxa(dat)[1]
    
    n.covered <- sum(dat$Tab[s,] > thres)
    present <- table(droplevels(dat$Map$Group[ dat$Tab[s,] > thres ]))
    
    cat("\t",s,"\t",n.covered,"\n")
    
    if(length(present) > 1 & all(present > 2)){
      res <- data.frame(A = x[1], B = x[2], Species = s)
    }else{
      res <- NULL
    }
    
    Res <- rbind(Res,res)
  }
  
  return(Res)
  
})

Res <- do.call(rbind, Res)
Res
write.table(Res,"comparisons.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


sort(table(Res$Species))
