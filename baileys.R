# what are Bailey's doing with the plot switch??
library(ggplot2)

source('data_functions.R')
theme_set(theme_bw())
cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")

dat = get_data(startdate = "2013-03-11",min_num_plots = 24,treatments = c('CC','EC','XC'))

# just PB abundance
pb = make_N_data(species = 'PB',dat)

# monthly mean
pbmean = aggregate(pb$n,by=list(censusdate = pb$censusdate, treatment=pb$treatment), FUN=mean)

# plot data
baileys = ggplot(pb, aes(x = censusdate, y =n, colour = treatment)) +
  geom_jitter(height=.1,width=.3) +
  geom_line(aes(x=censusdate,y=x,colour=treatment),data=pbmean,size=1) +  #geom_smooth(method='loess',se = TRUE) +
  theme(legend.position = 'right') +
  labs(y = 'Abundance', x = NULL) +
  geom_vline(xintercept=as.Date('2015-04-10')) +
  scale_colour_manual(name = 'Treatment', values = cbPalette,
                      breaks=c("CC","EC","XC"),
                      labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
  scale_fill_manual(name = 'Treatment', values = cbPalette, 
                    breaks=c("CC","EC","XC"),
                    labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) 
baileys
ggsave('BaileysAbundance.png',baileys,width=6,height=2)


# PB disappear from EC plots but not EE
# PF maybe a difference?
# PP, PE, RM, BA no difference
# PM, PL, RF, RO not enough data