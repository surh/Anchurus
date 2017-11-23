library(ggplot2)
library(reshape2)
library(GGally)

Tab <- read.table("Porphyromonas_sp_57899_mk_results.txt", sep = "\t", header = TRUE)
head(Tab)

# Number of genes
nrow(Tab)

# Plot histograms of all p-values
Pvals <- Tab[,c("gene","hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
Pvals <- melt(Pvals, id.vars = "gene", value.name = "p.value", variable.name = "Test")
p1 <- ggplot(Pvals, aes(x = p.value)) +
  facet_wrap(~ Test, nrow = 2) +
  geom_histogram(bins = 20) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
ggsave("all_pvals.hist.svg", p1, width = 8, height = 4)

# Plot only those with at least one fixed non-synonymous mutation
Pvals <- subset(Tab, Dn > 0)
nrow(Pvals)
Pvals <- Pvals[,c("gene","hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
Pvals <- melt(Pvals, id.vars = "gene", value.name = "p.value", variable.name = "Test")
p1 <- ggplot(Pvals, aes(x = p.value)) +
  facet_wrap(~ Test, nrow = 2) +
  geom_histogram(bins = 20) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
ggsave("withDn_pvals.hist.svg", p1, width = 8, height = 4)


### Barplot number of tests
Pvals <- Tab[,c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
Pvals <- apply(is.na(Pvals),2,function(x) table(factor(x,levels=c(TRUE,FALSE))))
Pvals <- melt(Pvals, value.name = "Count", varnames = c("is.na","test"))
p1 <- ggplot(Pvals, aes(x = test, y = Count)) + 
  geom_bar(aes(fill = is.na), stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 90))
p1
ggsave("sucess_tests.barplot.svg", p1, width = 6, height = 4)

### Barplot significant
Pvals <- Tab[,c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
Pvals <- data.frame(melt(colSums(Pvals < 0.05,na.rm = TRUE ),value.name = "Significant"))
Pvals$Test <- row.names(Pvals)
p1 <- ggplot(Pvals, aes(x = Test, y = Significant)) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))
p1
ggsave("significant_tests.barplot.svg", p1, width = 6, height = 4)

### Plot correlation all
#Pvals <- Tab[,c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
pval_columns <- c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")
pval_columns <- which(colnames(Tab) %in% pval_columns )
p1 <- ggpairs(subset(Tab, !is.na(g_none_p)),columns = pval_columns) + theme(axis.text.x = element_text(angle = 90))
p1
ggsave("correlation_pvals_all.svg", p1, width = 10, height = 10)

### Plot correlation with Dn
Pvals <- subset(Tab, Dn > 0)
#Pvals <- Pvals[,c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")]
pval_columns <- c("hg_p","hg_p_pseudo","g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo","g_yates_p_pseudo","g_williams_p_pseudo")
pval_columns <- which(colnames(Pvals) %in% pval_columns )
p1 <- ggpairs(subset(Pvals, !is.na(g_none_p)),columns = pval_columns) + theme(axis.text.x = element_text(angle = 90))
p1
ggsave("correlation_pvals_withDn.svg", p1, width = 10, height = 10)

## Plot neutrality index vs P-value
# p1 <- ggplot(Tab, aes(x = log10(ratio_pseudo), y = -log10(g_williams_p_pseudo))) +
#   geom_hline(yintercept = -log10(0.05)) +
#   geom_point(aes(col = Dn + 1)) +
#   scale_color_gradient2(trans = "log2",low = "#e08214",mid = "#f7f7f7",high = "#8073ac", midpoint = log2(2)) +
#   theme(panel.background = element_blank(),
#         panel.grid = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA))
# p1
Dat <- Tab[,c("gene","ratio_pseudo","Dn","hg_p","hg_p_pseudo",
              "g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo",
              "g_yates_p_pseudo","g_williams_p_pseudo")]
Dat <- melt(Dat, id.vars = c("gene","ratio_pseudo","Dn"), value.name = "p.value", variable.name = "Test")
head(Dat)
p1 <- ggplot(Dat, aes(x = log10(ratio_pseudo), y = -log10(p.value))) +
  facet_wrap(~ Test, nrow = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_point(aes(col = Dn + 1)) +
  scale_color_gradient2(trans = "log2",low = "#e08214",mid = "#f7f7f7",high = "#8073ac", midpoint = log2(2)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
ggsave("all_pvals.ratio.svg", p1, width = 8, height = 4)

## Plot neutrality index vs P-value for Dn > 0
Dat <- Tab[,c("gene","ratio_pseudo","Dn","hg_p","hg_p_pseudo",
              "g_none_p","g_yates_p","g_williams_p","g_none_p_pseudo",
              "g_yates_p_pseudo","g_williams_p_pseudo")]
Dat <- melt(Dat, id.vars = c("gene","ratio_pseudo","Dn"), value.name = "p.value", variable.name = "Test")
Dat <- subset(Dat, Dn > 0)
head(Dat)
p1 <- ggplot(Dat, aes(x = log10(ratio_pseudo), y = -log10(p.value))) +
  facet_wrap(~ Test, nrow = 2) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_point(aes(col = Dn + 1)) +
  scale_color_gradient2(trans = "log2",low = "#e08214",mid = "#f7f7f7",high = "#8073ac", midpoint = log2(3)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
p1
ggsave("withDn_pvals.ratio.svg", p1, width = 8, height = 4)

## Plot per position
Pvals <- Tab
Pvals$Pos <- (Pvals$start + Pvals$end) / 2
Pvals <- melt(Pvals,id.vars = c("gene","contig","start","end","Dn","Ds","Pn","Ps",
                                "ni","ratio","ratio_pseudo","hg_odds", "hg_odds_pseudo",
                                "alpha","alpha_pseudo","Pos"), variable.name = "Test",
              value.name = "p.value")
p1 <- ggplot(Pvals,aes(x = Pos, y = ni)) +
  facet_grid(Test ~ contig, scales = "free_x", space = "free_x") +
  geom_point(aes(size = Dn + 1, col = p.value < 0.05)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 360))
p1
ggsave("all_mk_by_position.svg",p1, width = 15, height = 7)

dat <- subset(Pvals, Test == "hg_p")
dat <- droplevels(dat)
p1 <- ggplot(dat,aes(x = Pos, y = log10(ratio_pseudo))) +
  facet_grid( ~ contig, scales = "free_x", space = "free_x") +
  geom_point(aes(size = Dn + 1, col = p.value < 0.05)) +
  scale_size_continuous(labels = c(0,1,2,4,8,12), breaks = c(1,2,3,5,9,13),name = "Dn") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_text(face = "bold", size = 20, color = "black"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))
p1
ggsave("all_mk_hg_position.svg",p1, width = 18, height = 5)
ggsave("all_mk_hg_position.png",p1, width = 18, height = 5)

Pvals <- subset(Pvals, Dn > 0)
p1 <- ggplot(Pvals,aes(x = Pos, y = ni)) +
  facet_grid(Test ~ contig, scales = "free_x", space = "free_x") +
  geom_point(aes(size = Dn + 1, col = p.value < 0.05)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 360))
p1
ggsave("withDn_mk_by_position.svg",p1, width = 15, height = 7)

### Get annotations
Genes <- read.table("genome.features", sep = "\t", header = TRUE)
head(Genes)
Pvals <- Tab
Pvals$n.significant <- rowSums(Pvals[,c("hg_p","hg_p_pseudo",
                                        "g_none_p","g_yates_p",
                                        "g_williams_p","g_none_p_pseudo",
                                        "g_yates_p_pseudo","g_williams_p_pseudo")] < 0.05,na.rm = TRUE)
# Pvals <- droplevels(subset(Pvals, n.significant > 0))
# Pvals <- Pvals[ order(Pvals$contig, Pvals$start), ]
Pvals <- droplevels(subset(Pvals, hg_p < 0.05))
Pvals <- Pvals[ order(Pvals$hg_p, decreasing = FALSE), ]
head(Pvals)

row.names(Genes) <- as.character(Genes$gene_id)
Pvals$functions <- as.character(Genes[as.character(Pvals$gene), "functions"])

Pvals
head(Pvals)
Pvals <- Pvals[,c("gene","contig","start","end","Dn","Ds","Pn","Ps","hg_p","n.significant","functions")]
write.table(Pvals,"significant.txt", col.names = TRUE, sep = "\t", quote = FALSE, row.names = FALSE)
