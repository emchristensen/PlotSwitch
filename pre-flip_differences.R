# ========================================================================
# kruskal-wallis test ----
# EMC 8/13/18
# we just want to know if there are differences in treatments before flip. we don't need to model trends/seasonality
# two sample t-test wont work -- data is not normal
# Used Kruskal-Wallis to test for significant difference, and pairwise wilcox to see what the differences are
# http://www.sthda.com/english/wiki/kruskal-wallis-test-in-r
# I also used k-w to test for differences between plots within treatments. none were (removals were almost significant for energy and dipos)

#' @description prepare data for doing rmanova
#' @param df data frame containing censusdate and plot as columns
#' 
prepare_dataframe = function(df) {
  df$censusdate = as.Date(df$censusdate)
  timecol = data.frame(censusdate = unique(df$censusdate))
  timecol$ordereddate = as.numeric(row.names(timecol))
  dat = merge(df,timecol)
  df$plot = as.factor(df$plot)
  return(dat)
}

cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")

# dipos ---- 
df1 = read.csv('Dipo_counts.csv')
dat1 = prepare_dataframe(df1) %>% filter(ordereddate<22)

# check that data within factor levels is normal
dip_control = dplyr::filter(dat1,treatment=='CC')
dip_excl = dplyr::filter(dat1,treatment=='EC')
dip_rem = dplyr::filter(dat1,treatment=='XC')

shapiro.test(dip_control$n)
# none of the treatments is normal. need a kruskal-wallis test 
kruskal.test(n~treatment, data=dat1)
# ok there are differences. use pairwise wilcox to see what differences are
pairwise.wilcox.test(dat1$n, dat1$treatment,
                     p.adjust.method = "bonferroni")
# there are pairwise differences between all 3 treatments

dipobox = ggplot(dat1,aes(x=treatment,y=n, colour=treatment)) + 
  geom_boxplot() +
  ggtitle('Kangaroo rat abundances 2013-2015') +
  scale_colour_manual(values=cbPalette) +
  labs(y = 'Count', x = NULL) +
  scale_x_discrete(breaks=c("CC","EC","XC"),
                   labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
  theme(legend.position = 'none')
dipobox

dat1 %>%
  group_by(treatment) %>%
  summarise(n = mean(n), sd = sd(n))

ggsave('DipoCountBoxplot.png',dipobox,width=6,height=2.5)

# small granivores ----
df = read.csv('SmallGranivores.csv')
dat = prepare_dataframe(df) %>% filter(ordereddate<22)

# check that data within factor levels is normal
sm_control = dplyr::filter(dat,treatment=='CC')
sm_excl = dplyr::filter(dat,treatment=='EC')
sm_rem = dplyr::filter(dat,treatment=='XC')

shapiro.test(sm_rem$n)
# exclosures are normal, removals and controls are not. need a kruskal-wallis test 
kruskal.test(n~treatment, data=dat)
# ok there are differences. use pairwise wilcox to see what differences are
pairwise.wilcox.test(dat$n, dat$treatment,
                     p.adjust.method = "bonferroni")
# there are pairwise differences between all 3 treatments

smgranbox = ggplot(dat,aes(x=treatment,y=n, colour=treatment)) + 
  geom_boxplot() +
  ggtitle('Small granivore abundances 2013-2015') +
  scale_colour_manual(values=cbPalette) +
  labs(y = 'Count', x = NULL) +
  scale_x_discrete(breaks=c("CC","EC","XC"),
                   labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
  theme(legend.position = 'none')
smgranbox

dat %>%
  group_by(treatment) %>%
  summarise(n = mean(n), sd = sd(n))

ggsave('SmallGranivoreBoxplot.png',smgranbox,width=6,height=2.5)

# energy ---- 
df2 = read.csv('TotalCommunityEnergy.csv')
dat2 = prepare_dataframe(df2) %>% filter(ordereddate<22)

# check that data within factor levels is normal
en_control = dplyr::filter(dat2,treatment=='CC')
en_excl = dplyr::filter(dat2,treatment=='EC')
en_rem = dplyr::filter(dat2,treatment=='XC')

shapiro.test(en_rem$n)
# exclosures and control are normal, removals are not. need a kruskal-wallis test 
kruskal.test(n~treatment, data=dat2)
# ok there are differences. use pairwise wilcox to see what differences are
pairwise.wilcox.test(dat2$n, dat2$treatment,
                     p.adjust.method = "bonferroni")
# there are pairwise differences between all 3 treatments

energybox = ggplot(dat2,aes(x=treatment,y=n, colour=treatment)) + 
  geom_boxplot() +
  ggtitle('Metabolic Flux 2013-2015') +
  scale_colour_manual(values=cbPalette) +
  labs(y = 'Metabolic Flux', x = NULL) +
  scale_x_discrete(breaks=c("CC","EC","XC"),
                   labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
  theme(legend.position = 'none')
energybox

dat2 %>%
  group_by(treatment) %>%
  summarise(n = mean(n), sd = sd(n))

ggsave('TotalEnergyBoxplot.png',energybox,width=6,height=2.5)




# differences between plots within treatment ----
# dipos
kruskal.test(n~plot,data=dip_control)
kruskal.test(n~plot,data=dip_excl)
kruskal.test(n~plot,data=dip_rem)

# small granivores
kruskal.test(n~plot,data=sm_control)
kruskal.test(n~plot,data=sm_excl)
kruskal.test(n~plot,data=sm_rem)

# metabolic flux
kruskal.test(n~plot,data=en_control)
kruskal.test(n~plot,data=en_excl)
kruskal.test(n~plot,data=en_rem)


# species richness ---- 
# df3 = read.csv('SpeciesRichness.csv')
# dat3 = prepare_dataframe(df3) %>% filter(ordereddate<22)
# 
# # check that data within factor levels is normal
# sr_control = dplyr::filter(dat3,treatment=='CC')
# sr_excl = dplyr::filter(dat3,treatment=='EC')
# sr_rem = dplyr::filter(dat3,treatment=='XC')
# 
# shapiro.test(sr_rem$n)
# # exclosures and control are normal, removals are not. need a kruskal-wallis test 
# kruskal.test(n~treatment, data=dat3)
# # ok there are differences. use pairwise wilcox to see what differences are
# pairwise.wilcox.test(dat3$n, dat3$treatment,
#                      p.adjust.method = "bonferroni")
# # there are pairwise differences between all 3 treatments
# 
# richbox = ggplot(dat3,aes(x=treatment,y=n, colour=treatment)) + 
#   geom_boxplot() +
#   ggtitle('Species Richness 2013-2015') +
#   scale_colour_manual(values=cbPalette) +
#   labs(y = '# Species', x = NULL) +
#   scale_x_discrete(breaks=c("CC","EC","XC"),
#                    labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
#   theme(legend.position = 'none')
# richbox
# 
# dat3 %>%
#   group_by(treatment) %>%
#   summarise(n = mean(n), sd = sd(n))
# 
# ggsave('SpeciesRichnessBoxplot.png',richbox,width=6,height=2.5)
