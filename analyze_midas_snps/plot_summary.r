library(ggplot2)

dat <- read.table("summary.txt")
head(dat)
colnames(dat) <- c("Strain", "Comparison", 'tested_genes', 'genes_with_Dn', 'pval0.1', 'pval0.05')
head(dat)


dat$E0.1 <- dat$tested_genes*0.1
dat$E0.05 <- dat$tested_genes*0.05

head(dat)

dat2 <- reshape2::melt(dat,id.vars = c("Strain", "Comparison", "tested_genes","genes_with_Dn"))
head(dat2)
dat0.1 <- droplevels(subset(dat2, variable %in% c("pval0.1","E0.1")))
dat0.1$variable <- as.numeric(dat0.1$variable)
dat0.05 <- droplevels(subset(dat2, variable %in% c("pval0.05","E0.05")))
dat0.05$variable <- as.numeric(dat0.05$variable)

p1 <- ggplot(dat2) +
  geom_line(data = dat0.1,aes(x = variable, y = value, group = Comparison:Strain), color = "red") +
  geom_line(data = dat0.05,aes(x = variable, y = value, group = Comparison:Strain), color = "blue") +
  theme_classic()
p1
ggsave("significant_vs_expected.svg",p1, width = 6, height = 4)
