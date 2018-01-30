library(ggplot2)

dir <- opts[1]
dir <- "Porphyromonas_sp_57899/"
file <- "Porphyromonas_sp_57899/mk_results.Supragingival.plaque_Tongue.dorsum.txt"
file <- "Streptococcus_infantis_62471/mk_results.Buccal.mucosa_Tongue.dorsum.txt"


Dat <- read.table(file,header = TRUE, sep = "\t")
head(Dat)

summary(Dat$Dn)
p1 <- ggplot(Dat, aes(x = Dn)) +
  geom_histogram(bins = 9) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_log10(breaks = c(1, 10, 50, 150, 1000)) +
  AMOR::theme_blackbox
p1
ggsave("nD_histogram.svg",p1, width = 6, height = 4)

p1 <- ggplot(Dat, aes(x = hg_p)) +
  geom_histogram(bins = 20) +
  AMOR::theme_blackbox
p1
ggsave("hgP_histogram.svg",p1, width = 6, height = 4)

subset(Dat, hg_p < 0.05)

## Plot per position
Pvals <- Dat
Pvals$Pos <- (Pvals$start + Pvals$end) / 2
Pvals <- reshape2::melt(Pvals,id.vars = c("gene","contig","start","end","Dn","Ds","Pn","Ps",
                                "ni","ratio","ratio_pseudo","hg_odds", "hg_odds_pseudo",
                                "alpha","alpha_pseudo","Pos"), variable.name = "Test",
              value.name = "p.value")
Pvals <- droplevels(subset(Pvals, Test == 'hg_p'))
head(Dat)
p1 <- ggplot(Pvals,aes(x = Pos, y = ni)) +
  facet_grid(. ~ contig, scales = "free_x", space = "free_x") +
  geom_point(aes(size = Dn + 1, col = p.value < 0.05)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 360))
p1
ggsave("mk_by_position.svg",p1, width = 10, height = 3)

p1 <- ggplot(Pvals,aes(x = Pos, y = -log10(p.value) * sign(ni))) +
  facet_grid(. ~ contig, scales = "free_x", space = "free_x") +
  geom_point(aes(size = Dn + 1, col = p.value < 0.05)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle = 360))
p1
ggsave("mk_by_position_logp.svg",p1, width = 10, height = 3)
