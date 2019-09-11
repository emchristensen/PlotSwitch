# comparing values (number of dipos, number of non-dipo species, and total energy) between treatments, before treatment change in 2015

library(multcomp)
library(dplyr)
library('mgcv')
library(ggplot2)

cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")
# dipos -----
df1 = read.csv('Data/Dipo_counts.csv') %>% 
  filter(numericdate<16.543) %>%
  mutate(plot = factor(plot))
# plot data
ggplot(df1, aes(x=censusdate, y=n, colour=treatment)) + geom_jitter(width=.1, height=.2)

# model krat abundance
dipo.mod <- gam(n ~ treatment + numericdate + s(plot, bs = 're'), data = df1,
                family = poisson(), method ='REML')
summary(dipo.mod)

## Achim Zeileis answered a Q on Crossvalidated on this topic, so we can
## follow him and use a hand-crafted contrast matrix.
## https://stats.stackexchange.com/a/376257/1390

# test for pairwise differences using glht (general linear hypothesis)
contr <- matrix(0, nrow = 3, ncol = length(coef(dipo.mod)))
colnames(contr) <- names(coef(dipo.mod))
rownames(contr) <- c("EC - CC", "XC - CC", "XC - EC")
contr[, 2:3] <- rbind(c(1, 0), c(0, 1), c(-1, 1))
dipo.glht <- glht(dipo.mod, linfct = contr)
summary(dipo.glht)                        # summary

# small granivores ----
df2 = read.csv('Data/SmallGranivores.csv') %>% 
  filter(numericdate<16.543) %>%
  mutate(plot = factor(plot))
ggplot(df2, aes(x=censusdate, y=n, colour=treatment)) + geom_jitter(width=.1, height=.2)

# model small granivore abundance
sm.mod <- gam(n ~ treatment + numericdate + s(plot, bs = 're'), data = df2,
              family = poisson(), method = 'REML')
summary(sm.mod)
plot(sm.mod, pages = 1)

# test for pairwise differences using glht (general linear hypothesis)
contr <- matrix(0, nrow = 3, ncol = length(coef(sm.mod)))
colnames(contr) <- names(coef(sm.mod))
rownames(contr) <- c("EC - CC", "XC - CC", "XC - EC")
contr[, 2:3] <- rbind(c(1, 0), c(0, 1), c(-1, 1))
sm.glht <- glht(sm.mod, linfct = contr)
summary(sm.glht)                        # summary

# energy ---- 
df3 = read.csv('Data/TotalCommunityEnergy.csv') %>% 
  filter(numericdate<16.543) %>%
  mutate(plot = factor(plot))
ggplot(df3, aes(x=censusdate, y=n, colour=treatment)) +geom_jitter(width=.1, height=.2)

# model total metabolic energy
en.mod <- gam(n ~ treatment + numericdate + s(plot, bs = 're'), data = df3,
              family = tw, method = "REML")
summary(en.mod)

# test for pairwise differences using glht (general linear hypothesis)
contr <- matrix(0, nrow = 3, ncol = length(coef(en.mod)))
colnames(contr) <- names(coef(en.mod))
rownames(contr) <- c("EC - CC", "XC - CC", "XC - EC")
contr[, 2:3] <- rbind(c(1, 0), c(0, 1), c(-1, 1))
en.glht <- glht(en.mod, linfct = contr)
summary(en.glht)                        # summary


# ========================================================================
# create boxplots for supplement
dipobox = ggplot(df1,aes(x=treatment,y=n, colour=treatment)) + 
  geom_boxplot() +
  ggtitle('Kangaroo rat abundances 2013-2015') +
  scale_colour_manual(values=cbbPalette[c(6,1,4)]) +
  labs(y = 'Count', x = NULL) +
  scale_x_discrete(breaks=c("CC","EC","XC"),
                   labels=c("Control", "Kangaroo rat+", "Rodent+")) +
  theme(legend.position = 'none')
dipobox
#ggsave('Figures/DipoCountBoxplot.pdf',dipobox,width=6,height=2.5)
#ggsave('Figures/DipoCountBoxplot.tiff',dipobox,width=6,height=2.5)

smgranbox = ggplot(df2,aes(x=treatment,y=n, colour=treatment)) + 
  geom_boxplot() +
  ggtitle('Small granivore abundances 2013-2015') +
  scale_colour_manual(values=cbbPalette[c(6,1,4)]) +
  labs(y = 'Count', x = NULL) +
  scale_x_discrete(breaks=c("CC","EC","XC"),
                   labels=c("Control", "Kangaroo rat+", "Rodent+")) +
  theme(legend.position = 'none')
smgranbox

df2 %>%
  group_by(treatment) %>%
  summarise(avg = mean(n), sd = sd(n))

ggsave('Figures/SmallGranivoreBoxplot.pdf',smgranbox,width=6,height=2.5)
ggsave('Figures/SmallGranivoreBoxplot.tiff',smgranbox,width=6,height=2.5)

energybox = ggplot(df3,aes(x=treatment,y=n, colour=treatment)) + 
  geom_boxplot() +
  ggtitle('Metabolic Flux 2013-2015') +
  scale_colour_manual(values=cbbPalette[c(6,1,4)]) +
  labs(y = 'Metabolic Flux', x = NULL) +
  scale_x_discrete(breaks=c("CC","EC","XC"),
                   labels=c("Control", "Kangaroo rat+", "Rodent+")) +
  theme(legend.position = 'none')
energybox

df3 %>%
  group_by(treatment) %>%
  summarise(avg = mean(n), sd = sd(n))

#ggsave('Figures/TotalEnergyBoxplot.pdf',energybox,width=6,height=2.5)
#ggsave('Figures/TotalEnergyBoxplot.tiff',energybox,width=6,height=2.5)
