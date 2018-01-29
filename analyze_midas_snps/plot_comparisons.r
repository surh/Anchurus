library(ggplot2)

dat <- read.table("comparisons.txt", sep = "\t", header = TRUE)
head(dat)
dat$Comparison <- paste(dat$A, dat$B, sep = ".")
dat$Value <- 1
dat$Species <- factor(dat$Species, levels = names(sort(table(dat$Species))))
head(dat)


p1 <- ggplot(dat,aes(x = Species, y = Comparison)) +
  geom_tile(aes(fill = Value)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
p1

ggsave("comparsions.svg", p1, width = 12, height = 4)



