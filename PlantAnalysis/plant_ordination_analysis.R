# code for doing ordination analysis on Portal plant communities (summer and winter annuals)
# I followed Sarah's code from Supp et al 2012 https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/12-0370.1
# Sarah did a pcca on the plant community (summer and winter annuals)

# another example: https://rgriff23.github.io/2017/05/23/mosquito-community-ecology-in-vegan.html



library(dplyr)
library(vegan)
library(ggplot2)
library(cowplot)

cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")
theme_set(theme_bw())

# functions ----
#' @title plant pcca
#' 
#' @description prepares plant data for pcca and returns cca object
#' @param plantdat matrix of species counts by plot and year
#'
plant_pcca = function(plantdat) {
  # sqrt transofrm abundance data to account for huge differences in abundance year to year
  species = sqrt(as.matrix(plantdat[,!names(plantdat) %in% c('year','season','plot','flip_type','treat_before','treat_after')]))
  year = as.factor(plantdat$year)
  trt = as.factor(plantdat$treat_before)
  plot = as.factor(plantdat$plot)
  
  # pcca (constrained by year)
  plant.pcca <- cca(species ~ trt + Condition(year))
  
  return(plant.pcca)
}

#' @title plot pcca with ellipses
#' @param pcca.obj result of pcca analysis
#' @param plantdat original matrix of plant count data used in pcca analysis
#' @param title string; title of plot
#' @param Palette vector of hex colors
#' 
plot_pcca_ellipses = function(pcca.obj, plantdat, title, Palette) {
  
  results=data.frame(scores(pcca.obj)$sites,
                      year=plantdat$year,
                      plot=plantdat$plot,
                      treatment=plantdat$treat_before)
  plotobj=ggplot(results,aes(CCA1, CA1)) + 
    stat_ellipse(aes(color = treatment)) +
    geom_point(aes(color = treatment)) +
    #geom_text(aes(label = plot, color = as.factor(year))) +
    coord_equal() +                     # biplots only work with equal scaling
    scale_shape_discrete(guide=F) +
    theme(legend.title=element_blank()) +
    #scale_color_discrete(direction=-1) +
    scale_colour_manual(values = Palette,
                        breaks=c("control","exclosure","removal"),
                        labels=c("long-term\ncontrol\n", "kangaroo rat\nremoval\n", "rodent\nremoval")) +
    ggtitle(title)
  return(plotobj)
}

# data ----
# Data I'm using is all censuses before flip, back to 2008 (many summers were missed)
# Summer: 2008, 2011, 2014 and Winter: 2008, 2012, 2013, 2014, 2015

dat.winter <- read.csv('WinterAnnualTreatments.csv')
dat.summer <- read.csv('SummerAnnualTreatments.csv')

# ===================================================================
# Do plant communities differ by treatment 2008-2015?

# krat exclosures vs removals ----
dat.winter1 = dplyr::filter(dat.winter,treat_before %in% c('exclosure','removal'))
dat.summer1 = dplyr::filter(dat.summer,treat_before %in% c('exclosure','removal'))

# pcca model: winter
win.pcca1 = plant_pcca(dat.winter1)

vif.cca(win.pcca1)
anova(win.pcca1)
permutest(win.pcca1,permutations=500) # should be similar to anova on pcca
#anova(win.pcca1,strata=dat.winter1$year) # more conservative test   

# proportion of variance explained
win.pcca1$CCA$tot.chi/win.pcca1$tot.chi

# plot with year, plot, treatment
excl_rem_win = plot_pcca_ellipses(win.pcca1, dat.winter1,'Winter Annual Plants',Palette=cbPalette[2:3])
excl_rem_win
#ggsave('Winter_Exclosure_Removal.png',excl_rem_win,width=4,height=3)

# pcca model: summer
sum.pcca1 = plant_pcca(dat.summer1)

vif.cca(sum.pcca1)
anova(sum.pcca1)
permutest(sum.pcca1,permutations=500) 

# proportion of variance explained
sum.pcca1$CCA$tot.chi/sum.pcca1$tot.chi

# plot with year, plot, treatment
excl_rem_sum = plot_pcca_ellipses(sum.pcca1, dat.summer1,'Summer Annual Plants',Palette=cbPalette[2:3])
excl_rem_sum
#ggsave('Summer_Exclosure_Removal.png',excl_rem_sum,width=4,height=3)

# cowplot grid
excl_rem_row <- plot_grid( excl_rem_win + theme(legend.position="none"),
                   excl_rem_sum + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("A", "B"),
                   hjust = -1,
                   nrow = 1)
legend1 <- get_legend(excl_rem_win)
excl_rem <- plot_grid( excl_rem_row, legend1, rel_widths = c(3, .6))
excl_rem
ggsave('Plants_Exclosure_Removal.png',excl_rem,width=8,height=3)


# controls vs total rodent removals ----
dat.winter3 = dplyr::filter(dat.winter,treat_before %in% c('removal','control'))
dat.summer3 = dplyr::filter(dat.summer,treat_before %in% c('removal','control'))

# pcca model: winter
win.pcca3 = plant_pcca(dat.winter3)

vif.cca(win.pcca3)
anova(win.pcca3)
permutest(win.pcca3,permutations=500) # should be similar to anova on pcca
#anova(win.pcca3,strata=dat.winter3$year) # more conservative test   

# proportion of variance explained
win.pcca3$CCA$tot.chi/win.pcca3$tot.chi

# plot with year, plot, treatment
ctrl_rem_win = plot_pcca_ellipses(win.pcca3, dat.winter3,'Winter Annual Plants',Palette=cbPalette[c(1,3)])
ctrl_rem_win
#ggsave('Winter_Control_Removal.png',ctrl_rem_win,width=4,height=3)

# pcca: summer
sum.pcca3 = plant_pcca(dat.summer3)

vif.cca(sum.pcca3)
anova(sum.pcca3)
permutest(sum.pcca3,permutations=500)
#anova(sum.pcca3,strata=dat.summer3$year) # more conservative test 

# proportion of variance explained
sum.pcca3$CCA$tot.chi/sum.pcca3$tot.chi

# plot with year, plot, treatment
ctrl_rem_sum = plot_pcca_ellipses(sum.pcca3, dat.summer3,'Summer Annual Plants',Palette=cbPalette[c(1,3)])
ctrl_rem_sum
#ggsave('Summer_Control_Removal.png',ctrl_rem_sum,width=4,height=3)

# cowplot grid
ctrl_rem_row <- plot_grid( ctrl_rem_win + theme(legend.position="none"),
                           ctrl_rem_sum + theme(legend.position="none"),
                           align = 'vh',
                           labels = c("A", "B"),
                           hjust = -1,
                           nrow = 1)
legend2 <- get_legend(ctrl_rem_win)
ctrl_rem <- plot_grid( ctrl_rem_row, legend2, rel_widths = c(3, .6))
ctrl_rem
ggsave('Plants_Control_Removal.png',ctrl_rem,width=8,height=3)


# controls vs krat exclosures ----
dat.winter2 = dplyr::filter(dat.winter,treat_before %in% c('exclosure','control'))
dat.summer2 = dplyr::filter(dat.summer,treat_before %in% c('exclosure','control'))

# pcca model: winter
win.pcca2 = plant_pcca(dat.winter2)

vif.cca(win.pcca2)
anova(win.pcca2)
permutest(win.pcca2,permutations=500) # should be similar to anova on pcca
#anova(win.pcca2,strata=dat.winter2$year) # more conservative test   

# proportion of variance explained
win.pcca2$CCA$tot.chi/win.pcca2$tot.chi

# plot with year, plot, treatment
ctrl_excl_win = plot_pcca_ellipses(win.pcca2, dat.winter2,'Winter Annual Plants',Palette=cbPalette[1:2])
ctrl_excl_win
#ggsave('Winter_Control_Exclosure.png',ctrl_excl_win,width=4,height=3)

# pcca: summer
sum.pcca2 = plant_pcca(dat.summer2)

vif.cca(sum.pcca2)
anova(sum.pcca2)
permutest(sum.pcca2,permutations=500)
#anova(sum.pcca2,strata=dat.summer2$year) # more conservative test 

# proportion of variance explained
sum.pcca2$CCA$tot.chi/sum.pcca2$tot.chi

# plot with year, plot, treatment
ctrl_excl_sum = plot_pcca_ellipses(sum.pcca2, dat.summer2,'Summer Annual Plants',Palette=cbPalette[1:2])
ctrl_excl_sum
#ggsave('Summer_Control_Exclosure.png',ctrl_excl_sum,width=4,height=3)

# cowplot grid
ctrl_excl_row <- plot_grid( ctrl_excl_win + theme(legend.position="none"),
                           ctrl_excl_sum + theme(legend.position="none"),
                           align = 'vh',
                           labels = c("A", "B"),
                           hjust = -1,
                           nrow = 1)
legend3 <- get_legend(ctrl_excl_win)
ctrl_excl <- plot_grid( ctrl_excl_row, legend3, rel_widths = c(3, .6))
ctrl_excl
ggsave('Plants_Control_Exclosure.png',ctrl_excl,width=8,height=3)


# other significance tests ----
### ADONIS - another test to look for compositional differences, with similar results to above.(from Supp et al 2012)
# winter
#win.spp.canb = vegdist(win.species, method = "canb")
#win.canb = adonis(win.spp.canb ~ win.trt, permutation=1000)
#win.canb
