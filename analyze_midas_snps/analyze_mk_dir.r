library(ggplot2)
library(reshape2)

indir <- 'out/'

species_dirs <- list.dirs(indir,recursive = FALSE)

Res <- NULL
for (s in species_dirs){
  #s <- species_dirs[2]
  
  files <- list.files(s)
  files <- files[ grep("mk_results.",x = files)] 
  #files
  
  for (f in files){
    #f <- files[1]
    
    info <- sub("mk_results.",replacement = "",x = f)
    info <- sub(".txt",replacement = "",x = info)
    info <- strsplit(x = info,split = "_")[[1]]
    #info
    
    cat(s,"\t",info,"\n")
    
    infile <- paste(s,"/",f,sep = "")
    #infile
    mk <- read.table(infile, header = TRUE, sep = "\t")
    #head(mk)

    if(nrow(mk) == 0){
      cat("\tSkip\n")
      next
    }
    
    mk$Pos <- (mk$start + mk$end) / 2
    
    # Plot
    p1 <- ggplot(mk,aes(x = Pos, y = log10(ratio_pseudo))) +
      facet_grid( ~ contig, scales = "free_x", space = "free_x") +
      geom_point(aes(size = Dn + 1, col = hg_p < 0.05)) +
      scale_size_continuous(labels = c(0,1,2,4,8,12), breaks = c(1,2,3,5,9,13),name = "Dn") +
      theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            axis.text.x = element_blank(),
            strip.text.x = element_blank(),
            axis.text.y = element_text(face = "bold", size = 20, color = "black"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18))
    #p1
    filename <-paste(basename(s),"_",info[1],"_",info[2],"_mk.genome.svg", sep = "")
    ggsave(filename,p1, width = 20, height = 5)
    
    res <- subset(mk, hg_p < 0.05 )
    if(nrow(res) == 0){
      cat("\tNO SIGNIFICANT\n")
      next
    }
    
    res$ni <- res$ratio <- res$hg_odds <- res$alpha <- res$alpha_pseudo <- NULL
    res$g_none_p <- res$g_none_p_pseudo <- res$g_williams_p <- res$g_williams_p_pseudo <- NULL
    res$hg_odds_pseudo <- res$hg_p_pseudo <- res$g_yates_p <- res$g_yates_p_pseudo <- NULL
    res$Species <- basename(s)
    res$group1 <- info[1]
    res$group2 <- info[2]
    #head(res)
    
    Res <- rbind(Res,res)
  }
}

write.table(Res, "signifincant_all.txt",sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
