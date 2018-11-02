# code for doing a permanova to look for differences in winter/summer plant communities based on rodent treatments
# following example from:
#     https://rpubs.com/collnell/manova 

library(vegan)
library(dplyr)
library(ggplot2)

# plotting function ----
plot_betadisper = function(betadisper.obj, colorpalette = cbPalette, title) {
  df = as.data.frame(betadisper.obj$vectors)
  df$treatment = betadisper.obj$group
  
  centroids = data.frame(betadisper.obj$centroids[,1:2])
  centroids$treatment = row.names(centroids)
  
  plotobj = ggplot(df,aes(PCoA1, PCoA2)) + 
    stat_ellipse(aes(color = treatment),level=.95) +
    geom_point(aes(colour = treatment), data = centroids, size = 5, stroke = 1, shape = 3) +
    geom_point(aes(color = treatment)) +
    coord_equal() +                      # biplots only work with equal scaling
    scale_shape_discrete(guide=F) +
    theme(legend.title=element_blank(), legend.position = "right",
          legend.key.height = unit(1, "cm"), legend.spacing = unit(1, "cm"),
          legend.key.width = unit(1, "cm")) +
    scale_colour_manual(values = colorpalette,
                        breaks=c("control","exclosure","removal"),
                        labels=c("long-term\ncontrol", "kangaroo rat\nremoval", "rodent\nremoval")) +
    ggtitle(title)
}

# data ----
# Data: censuses before flip, back to 2008 (many summers were missed)
#    Summer: 2008, 2011, 2014
#    Winter: 2008, 2012, 2013, 2014, 2015

# select data only from plots that were used for the GAM analysis (10 plots)
gamplots = c(4,5,6,7,11,13,14,17,18,24)
#treatplots = c(5,6,7,24,13,18)
dat.winter <- read.csv('PlantAnalysis/WinterAnnualTreatments.csv', stringsAsFactors = F) %>% 
  filter(plot %in% gamplots)
dat.summer <- read.csv('PlantAnalysis/SummerAnnualTreatments.csv', stringsAsFactors = F)%>%
  filter(plot %in% gamplots)
winter.matrix <- dat.winter[,4:33]
summer.matrix <- dat.summer[,4:35]

# sqrt transform to minimize influence of most abundant species
winter.mat <- sqrt(winter.matrix)
summer.mat <- sqrt(summer.matrix)

# calculate ecological distance
winter.dist <- vegdist(winter.mat, method = 'bray')
summer.dist <- vegdist(summer.mat, method = 'bray')

# perMANOVA
winter.div <- adonis2(winter.dist~treat_before, data = dat.winter, permutations = 999, method = 'bray', strata = c('plot','year'))
winter.div
summer.div <- adonis2(summer.dist~treat_before, data = dat.summer, permutations = 999, method = 'bray', strata = c('plot','year'))
summer.div
dispersion.w <- betadisper(winter.dist, group=dat.winter$treat_before)
permutest(dispersion.w)
dispersion.s <- betadisper(summer.dist, group=dat.summer$treat_before)
permutest(dispersion.s)

# plot 

# base graphics
plot(dispersion.w, hull=F, ellipse=T, main = 'Winter Annual Plants')

# ggplot
cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")
plot.winter = plot_betadisper(dispersion.w, cbPalette, title='Winter Annual Plants')
plot.winter
plot.summer = plot_betadisper(dispersion.s, cbPalette, title='Summer Annual Plants')
plot.summer
