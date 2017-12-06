library(AMOR)


##### Read tables ####
orig <- read.table("../2017-08-07.abundances/merge.species/relative_abundance.txt", row.names = 1, header = TRUE)
xac <- read.table("xac/temp/merge.species/relative_abundance.txt", row.names = 1, header = TRUE)
xad <- read.table("xad/temp/merge.species/relative_abundance.txt", row.names = 1, header = TRUE)

# merge
all(row.names(orig) == row.names(xac))
all(row.names(xac) == row.names(xad))

all <- cbind(orig,xac,xad)

#### Read metadata, from 2.merge_all_metadata.r ####

# Load metadata from metagenome samples
file <- "~/micropopgen/data/HMP/2017-07-17.all_runs/all_runs.txt"
all_runs <- read.table(file,header = TRUE,sep = "\t")
head(all_runs)

file <- "~/micropopgen/data/HMP/HMASM.csv"
hmasm_samples <- read.table(file,sep = ",", quote = '"', header = TRUE)
head(hmasm_samples)

file <- "~/micropopgen/data/HMP/HMASM-690.csv"
hmasm_samples_690 <- read.table(file,header = TRUE,sep = ",",quote = '"')
head(hmasm_samples_690)

file <- "~/micropopgen/data/HMP/hmp_catalogue_wgs_all_metadata.csv"
catalogue <- read.csv(file)
head(catalogue)

# The first thing is to identify which samples are in the original publication,
# and which passed QC

# Compare body sites from hmasm sets
sites1 <- as.character(hmasm_samples_690$Body.Site)
names(sites1) <- as.character(hmasm_samples_690$SRS.ID)
sites2 <- as.character(hmasm_samples$Body.Site)
names(sites2) <- as.character(hmasm_samples$SRS.ID)

head(sites1)
head(sites2)

shared <- intersect(names(sites1),names(sites2))

if(!all(sites1[ shared ] == sites2[ shared ])){
  stop("ERROR: Not all samples match on body site")
}else{
  cat("All samples match\n")
}

hmasm_samples$HMASM.QC <- !is.na(sites1[ as.character(hmasm_samples$SRS.ID) ])
table(hmasm_samples$HMASM.QC)
head(hmasm_samples)
#row.names(hmasm_samples) <- as.character(hmasm_samples$SRS.ID)

# Add information about puiblication and QC into main table
ftable(HMP.Isolation.Body.Site ~ HMP.Isolation.Body.Subsite, catalogue)
row.names(hmasm_samples) <- as.character(hmasm_samples$SRS.ID)
catalogue$HMASM <- as.character(catalogue$Sequence.Read.Archive.ID) %in% row.names(hmasm_samples)
catalogue$HMASM.QC <- hmasm_samples[ as.character(catalogue$Sequence.Read.Archive.ID), "HMASM.QC" ]
table(catalogue$HMASM, useNA = "always")
table(catalogue$HMASM.QC, useNA = "always")

all(as.character(all_runs$secondary_sample_accession) %in% as.character(catalogue$Sequence.Read.Archive.ID))

# Set to stool all GI samples
catalogue$HMP.Isolation.Body.Subsite <- as.character(catalogue$HMP.Isolation.Body.Subsite)
catalogue$HMP.Isolation.Body.Subsite[ catalogue$HMP.Isolation.Body.Site == "gastrointestinal_tract" ] <- "Stool"

rm(all_runs,hmasm,hmasm_samples,hmasm_samples_690,file,shared,sites1,sites2)

##### Add set of samples ####
setdiff(colnames(all),catalogue$Sequence.Read.Archive.ID)
catalogue$Sample.Set <- NA
catalogue$Sample.Set[ catalogue$Sequence.Read.Archive.ID %in% colnames(orig) ] <- "Pilot"
catalogue$Sample.Set[ catalogue$Sequence.Read.Archive.ID %in% colnames(xac) ] <- "Round2"
catalogue$Sample.Set[ catalogue$Sequence.Read.Archive.ID %in% colnames(xad) ] <- "Round2"

rm(xac,xad,orig)

catalogue <- subset(catalogue, Sequence.Read.Archive.ID %in% colnames(all))
row.names(catalogue) <- as.character(catalogue$Sequence.Read.Archive.ID)
catalogue <- catalogue[ colnames(all), ]

Dat <- create_dataset(all, catalogue)
Dat <- clean(Dat, verbose = TRUE)
Dat


Dat.pca <- PCA(Dat, cor = FALSE)
summary(Dat.pca)
p1 <- plotgg(Dat.pca, col = "HMP.Isolation.Body.Subsite")
dat <- p1$data

# All
p1 <- ggplot(dat, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = HMP.Isolation.Body.Subsite), shape = 21, col = "black", size = 3) +
  xlab(label = "PC1 (22.75%)") +
  ylab(label = "PC2 (9.44%)") +
  theme_blackbox
p1
ggsave("PCA_all_RA.svg", p1, width = 8, height = 6)

# Pilot
p1 <- ggplot(subset(dat,Sample.Set == "Pilot"), aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = HMP.Isolation.Body.Subsite), shape = 21, col = "black", size = 3) +
  xlab(label = "PC1 (22.75%)") +
  ylab(label = "PC2 (9.44%)") +
  theme_blackbox
p1
ggsave("PCA_pilot_RA.svg", p1, width = 8, height = 6)

# QC
p1 <- ggplot(dat, aes(x = PC1, y = PC2)) +
  geom_point(data = subset(dat, !(!is.na(HMASM.QC) & HMASM.QC == FALSE)), aes(fill = HMASM.QC), shape = 21,
             col = "black", size = 3) +
  scale_fill_manual(values = c("lightblue","grey")) +
  geom_point(data = subset(dat, (!is.na(HMASM.QC) & HMASM.QC == FALSE)), shape = 21, col = "black", size = 3, fill = "red") +
  xlab(label = "PC1 (22.75%)") +
  ylab(label = "PC2 (9.44%)") +
  theme_blackbox
p1
ggsave("PCA_all_QC_RA.svg", p1, width = 8, height = 6)


