# Create some files for testing analysis scripts
setwd("~/micropopgen/exp/2017/today3/")

info <- read.table("Haemophilus_parainfluenzae_62356/snps_info.txt", sep = "\t", header = TRUE)
head(info)
freq <- read.table("Haemophilus_parainfluenzae_62356/snps_freq.txt", header = TRUE)
head(freq)
depth <- read.table("Haemophilus_parainfluenzae_62356/snps_depth.txt", header = TRUE)
head(depth)


new.samples <- c("sample3","sample4")

set.seed(3)
n_sites <- nrow(freq)
prop_new_sites <- 0.8
n_new_sites <-  round(n_sites * prop_new_sites)
n_samples_orig <- ncol(freq) - 1
i <- 1
for(s in new.samples){
  # s <- new.samples[1]
  
  new_col <- freq[ ,i + 1 ]
  index <- sample(1:nrow(freq),size = n_new_sites, replace = FALSE)
  new_col[index] <- round(runif(n = n_new_sites),4)
  
  #new_depth <- rbinom(n = n_sites,size = 1,prob = 0.3) * rpois(n = n_sites,lambda = mean(depth[,i + 1]))
  new_depth <- round(depth[,i+1] + rnorm(n = n_sites, mean = 0, sd = 6))
  new_depth[ new_depth < 0 ] <- 0
  
  freq[,s] <- new_col
  depth[,s] <- new_depth
  
  i <- i + 1
  if(i > n_samples_orig) i <- 1
}

head(freq)
head(depth)

dir.create("snps")
write.table(freq,file = "snps/snps_freq.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(info,file = "snps/snps_info.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(depth,file = "snps/snps_depth.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
