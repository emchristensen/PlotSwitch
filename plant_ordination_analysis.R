# code for doing ordination analysis on Portal plant communities (summer and winter annuals)
# I followed Sarah's code from Supp et al 2012 https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/12-0370.1
# Sarah did a pcca on the plant community (summer and winter annuals)

# another example: https://rgriff23.github.io/2017/05/23/mosquito-community-ecology-in-vegan.html


library(dplyr)
library(vegan)
library(ggplot2)

# data
dat.winter <- read.csv('WinterAnnualTreatments2.csv')
dat.summer <- read.csv('SummerAnnualTreatments.csv')

# ===================================================================
# Do plant communities differ by treatment 2008-2015?
# 

# pcca model
# sqrt transofrm abundance data to account for huge differences in abundance year to year
win.species = sqrt(as.matrix(dat.winter[,!names(dat.winter) %in% c('year','season','plot','flip_type','treat_before','treat_after')]))

win.year = as.factor(dat.winter$year)
win.trt = as.factor(dat.winter$treat_before)
win.plot = as.factor(dat.winter$plot)

# pcca (constrained by year)
win.pcca <- cca(win.species ~ win.trt + Condition(win.year))

vif.cca(win.pcca)
anova(win.pcca)
permutest(win.pcca,permutations=500) 

# proportion of variance explained
win.pcca$CCA$tot.chi/win.pcca$tot.chi


### ADONIS - another test to look for compositional differences, with similar results to above.(from Supp et al 2012)
# winter
win.spp.canb = vegdist(win.species, method = "canb")
win.canb = adonis(win.spp.canb ~ win.trt, permutation=1000)
win.canb


# From the anova and permutest it seems that treatment is significant, but explains a tiny amount of variance.

# plot with year, plot, treatment
results=data.frame(scores(win.pcca)$sites,
                   year=win.year,
                   plot=win.plot,
                   treatment=win.trt)
ggplot(results,aes(CCA1, CCA2)) + 
  stat_ellipse(aes(color = treatment)) +
  geom_text(aes(label = plot, color = as.factor(year))) +
  coord_cartesian(xlim = c(-5, 4), ylim = c(-7, 7)) +
  theme(legend.title=element_blank()) +
  scale_color_discrete(direction=-1)


# ==========================================================================================================
# Summer annuals
# pcca model
# sqrt transofrm abundance data to account for huge differences in abundance year to year
sum.species = sqrt(as.matrix(dat.summer[,!names(dat.summer) %in% c('year','season','plot','flip_type','treat_before','treat_after')]))

sum.year = as.factor(dat.summer$year)
sum.trt = as.factor(dat.summer$treat_before)
sum.plot = as.factor(dat.summer$plot)

# pcca (constrained by year)
sum.pcca <- cca(sum.species ~ sum.trt + Condition(sum.year))

vif.cca(sum.pcca)
anova(sum.pcca)
permutest(sum.pcca,permutations=500) 

# proportion of variance explained
sum.pcca$CCA$tot.chi/sum.pcca$tot.chi


### ADONIS - another test to look for compositional differences, with similar results to above.(from Supp et al 2012)
# winter
sum.spp.canb = vegdist(sum.species, method = "canb")
sum.canb = adonis(sum.spp.canb ~ sum.trt, permutation=1000)
sum.canb


# From the anova and permutest it seems that treatment is significant, but explains a tiny amount of variance.

# plot with year, plot, treatment
results2=data.frame(scores(sum.pcca)$sites,
                   year=sum.year,
                   plot=sum.plot,
                   treatment=sum.trt)
ggplot(results2,aes(CCA1, CCA2)) + 
  stat_ellipse(aes(color = treatment)) +
  geom_text(aes(label = plot, color = as.factor(year))) +
  coord_cartesian(xlim = c(-5, 4), ylim = c(-7, 7)) +
  theme(legend.title=element_blank()) +
  scale_color_discrete(direction=-1)

