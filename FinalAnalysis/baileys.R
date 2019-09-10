# Code to analyze C. Baileyi response to plot switch
library(ggplot2)
library("mgcv")

source('Data/data_functions.R')
theme_set(theme_bw())
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")

dat = get_data(startdate = "2013-03-11",min_num_plots = 24,treatments = c('CC','EC','XC'))

# just PB abundance
pb = make_N_data(species = 'PB',dat)



## Alternate version using model estimates
## add variables for modelling
pb <- mutate(pb,
             year = as.numeric(format(censusdate, format = "%Y")),
             doy = as.numeric(format(censusdate, format = "%j")),
             nmonth = as.numeric(format(censusdate, format = "%m")))

## ## push the boundary knot locations to past the endpoint of the data
## ## allows different values for Dec and Jan
## knots <- list(nmonth = c(0.5, 12.5))

m.pois <- gam(n ~ s(numericdate, k = 20, by = treatment),
              data = pb, family = poisson,
              method = 'REML')

gam.check(m.pois, rep = 100)

m.nb <- gam(n ~ s(numericdate, k = 20, by = treatment),
            data = pb, family = nb,
            method = 'REML')

gam.check(m.nb, rep = 100)

## as before the NegBIn is ~= poisson
## the check out is OK esp given the low mean, and
## randomised quantile residuals from https://gist.github.com/dill/d0e4549cbe6c903f232ff5bf456dd1f6
## aren't so bad - for the NegBin as better_check only works with tw and nb

## predict to get trend lines for each group

## new data, predict every 7 days; plot looks smooth enough
pdata <- with(pb, expand.grid(censusdate = seq(min(censusdate), max(censusdate),
                                               by = 7),
                              treatment = c('CC', 'EC', 'XC')))
pdata <- mutate(pdata,
                numericdate = as.numeric(censusdate) / 1000)
## predict
pred <- predict(m.pois, newdata = pdata, se.fit = TRUE)
ilink <- family(m.pois)$linkinv         # grab the inverse link
## add predictions and a creible interval
pdata <- mutate(pdata,
                fit     = pred$fit, se = pred$se.fit,
                predict = ilink(fit),
                upper   = ilink(fit + (2 * se)),
                lower   = ilink(fit - (2 * se)))

## plot
baileys <- ggplot(pb, aes(x = censusdate, y = n, colour = treatment)) +
  geom_ribbon(aes(x = censusdate, ymin = upper, ymax = lower, fill = treatment),
              data = pdata, alpha = 0.4, inherit.aes = FALSE) +
  geom_jitter(height = 0.1, width = 0.3, size=.5) +
  geom_line(aes(x = censusdate, y = predict, colour = treatment), data = pdata, size=1) +
  theme(legend.position = 'bottom',
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        axis.text = element_text(size=8)) +
  labs(y = 'Abundance', x = NULL) +
  geom_vline(xintercept=as.Date('2015-04-10')) +
  scale_colour_manual(name = 'Treatment:', values = cbbPalette[c(6,1,4)],
                      breaks=c("CC","EC","XC"),
                      labels=c("Control", "Kangaroo rat+", "Rodent+")) +
  scale_fill_manual(name = 'Treatment:', values = cbbPalette[c(6,1,4)], 
                    breaks=c("CC","EC","XC"),
                    labels=c("Control", "Kangaroo rat+", "Rodent+")) 
baileys

ggsave('Figures/BaileysAbundanceGAM.pdf', baileys, width=4, height = 2, dpi=300)
ggsave('Figures/BaileysAbundanceGAM.tiff', baileys, width=4, height = 2, dpi=300)


# ============================================================================================
# Old code: plot just monthly mean data without GAM 
## monthly mean
# pbmean = aggregate(pb$n,by=list(censusdate = pb$censusdate, treatment=pb$treatment), FUN=mean)
# 
# # plot data
# baileys = ggplot(pb, aes(x = censusdate, y =n, colour = treatment)) +
#   geom_jitter(height=.1,width=.3) +
#   geom_line(aes(x=censusdate,y=x,colour=treatment),data=pbmean,size=1) +  #geom_smooth(method='loess',se = TRUE) +
#   theme(legend.position = 'right') +
#   labs(y = 'Abundance', x = NULL) +
#   geom_vline(xintercept=as.Date('2015-04-10')) +
#   scale_colour_manual(name = 'Treatment', values = cbPalette,
#                       breaks=c("CC","EC","XC"),
#                       labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
#   scale_fill_manual(name = 'Treatment', values = cbPalette, 
#                     breaks=c("CC","EC","XC"),
#                     labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) 
# baileys
# ggsave('BaileysAbundance.png',baileys,width=6,height=2)
