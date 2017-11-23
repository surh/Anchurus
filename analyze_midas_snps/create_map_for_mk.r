
summary_file <- opts[1]
#summary_file <- "test/snps_summary.txt"

metadata_file <- "~/micropopgen/exp/2017/2017-07-17.metadata/samples_to_start.txt"

selected <- read.table(summary_file, header = TRUE)
meta <- read.table(metadata_file, header = TRUE, sep = "\t")
head(meta)

samples <- as.character(selected$sample_id)
meta <- subset(meta,Sequence.Read.Archive.ID %in% samples)

if(nrow(meta) != length(samples)){
  stop("ERROR")
}
  

meta <- data.frame(ID = meta$Sequence.Read.Archive.ID, Group = meta$HMP.Isolation.Body.Subsite)
table(meta$Group)
write.table(meta,"map.txt",sep = "\t", row.names = FALSE, quote = FALSE)
