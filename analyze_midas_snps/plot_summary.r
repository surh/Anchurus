library(ggplot2)

dat <- read.table("summary.txt", header = TRUE)
head(dat)
dat <- droplevels(subset(dat, ngenes > 0))
dat$Strain
nrow(dat)

dat$E0.1 <- dat$ngenes*0.1
dat$E0.05 <- dat$ngenes*0.05
head(dat)


dat2 <- dat[,c('Strain','Comparison','ngenes','n_Dn','n_Ds','n_D','n_P')]
dat2 <- reshape2::melt(dat2, id.vars = c('Strain','Comparison','ngenes'))
dat2$Percent.genes <- 100 * dat2$value / dat2$ngenes
head(dat2)

p1 <- ggplot(dat2, aes(x = variable, y = Percent.genes)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), shape = 19, size = 3) +
  theme_classic()
p1
ggsave("percent_mksites.svg", p1, width = 4, height = 4)

p1 <- ggplot(dat2, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_point(position = position_jitter(width = 0.2), shape = 19, size = 3) +
  theme_classic()
p1
ggsave("n_mksites.svg", p1, width = 4, height = 4)


p1 <- ggplot(dat, aes(x = 100 * n_Dn / ngenes, y = p0.1)) +
  geom_point(size = 3, color = 'red') +
  geom_smooth(method = "lm") +
  ylab("# p-value < 0.1") +
  theme_classic()
p1
ggsave("percentDN_vs_sig0.1.svg", p1, width = 4, height = 4)

p1 <- ggplot(dat, aes(x = 100 * n_Dn / ngenes, y = p0.05)) +
  geom_point(size = 3, color = 'blue') +
  geom_smooth(method = "lm") +
  ylab("# p-value < 0.05") +
  theme_classic()
p1
ggsave("percentDN_vs_sig0.05.svg", p1, width = 4, height = 4)


p1 <- ggplot(dat, aes(x = n_D, y = n_P)) +
  geom_point(size = 3, color = 'black') +
  geom_smooth(method = "lm") +
  theme_classic()
p1
ggsave("nD_vs_nP.svg", p1, width = 4, height = 4)


dat2 <- dat
dat2$E0.1 <- 0.1 * dat2$ngenes
dat2$E0.05 <- 0.05 * dat2$ngenes
head(dat2)

dat2 <- reshape2::melt(dat2,id.vars = c("Strain", "Comparison", "ngenes","n_Dn",'n_Ds','n_D','n_P'))
dat2$which <- "observed"
dat2$which[ dat2$variable %in% c('E0.1', 'E0.05')] <- 'expected'
dat2$which <- factor(dat2$which, levels = c('observed', 'expected'))
dat2$alpha <- 0.1
dat2$alpha[ dat2$variable %in% c('p0.05', 'E0.05')] <- 0.05
dat2$alpha <- factor(dat2$alpha)

head(dat2)
p1 <- ggplot(dat2, aes(x = which, y = value, group = Comparison:Strain:alpha, color = alpha)) +
  geom_line() +
  scale_color_manual(values = c('blue', 'red')) +
  theme_classic()
p1
ggsave("significant_vs_expected.svg",p1, width = 6, height = 4)
