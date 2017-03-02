# Questions: does # of Diops differ between CC, XC, and EC switched plots
#            does total energy of community differ?


library(vegan)
library(reshape)
library(ggplot2)
library(dplyr)
library(RCurl)


source('summarySE.r')
source('rodent_abundance_by_period_and_plot.r')


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

byspecies = rodent_abundance('Granivore')
treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                     'XC','EC','XC','CE',
                                     'CX','XX','CC','CX',
                                     'EC','CC','EE','XX',
                                     'CC','EC','EE','EE',
                                     'EE','CE','XX','XC'),plot=seq(1,24))
http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"
trapping = read.csv(text=getURL(http)) 
#trapping = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent_trapping.csv')
trapping$date = as.Date(paste(trapping$Year,trapping$Month,trapping$Day,sep='-'))
plotstrapped = aggregate(trapping$Sampled,by=list(period=trapping$Period),FUN=sum)
fullcensus = plotstrapped[plotstrapped$x>20,]

# ==========================================
byspecies$type = paste(byspecies$type_before,byspecies$type_after,sep='')

rdat = byspecies %>% filter(period>433,type %in% c('EC','CC','XC'),period %in% fullcensus$period) %>% select(plot,period,species,x,type)

# =============================================
# number of dipos per plot

d = filter(rdat,species %in% c('DM','DO','DS'))
dipos = aggregate(d$x,by=list(period=d$period,treatment=d$type,plot=d$plot),FUN=sum)
allplotsperiod = expand.grid(period=unique(dipos$period),plot=unique(dipos$plot))
allplotsperiod = merge(allplotsperiod,treatment)

dipos = merge(allplotsperiod,dipos,all=T)
dipos$x[is.na(dipos$x)] = 0

dipo_mn = aggregate(dipos$x,by=list(period=dipos$period,treatment=dipos$treatment),FUN=mean)
dpct <- summarySE(dipos, measurevar="x", groupvars=c("treatment","period"))


ggplot(dpct, aes(x=period, y=x, colour=treatment)) + 
  geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.1) +
  geom_line(size=1.5) +
  geom_point() +
  scale_y_continuous(name='# Dipos') +
  scale_x_continuous(name='') +
  geom_vline(xintercept = as.numeric(436),size=1.5) +
  scale_colour_manual(name="Experimental\nTreatment",
                      breaks=c("CC", "EC", "XC"),
                      labels=c("Long-term\nControl\n", "Established s.g.\ncommunity\n", "No established\ncommunity\n"),
                      values=cbPalette[c(6,4,7)])


# ======================================================================================
# total energy per plot
# my color scheme: CC = 6; EE=8; EC=4; CE=3; XC=7

allrodents = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','RO','SF','SH')
granivores  = c('BA','DM','DO','DS','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')
smgran = c('BA','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')

alldat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent.csv',stringsAsFactors = F)
source('rodent_energy_from_wgt.r')

rmass = alldat %>% filter(period>420,period %in% fullcensus$period, 
               species %in% smgran,
               !(plot %in% c(10,16,23)), note5 != 'R') %>%
  select(period,plot,species,wgt)

renergy = metabolic_energy_from_wgt(rmass)

energy = aggregate(renergy$energy,by=list(period=renergy$period,plot=renergy$plot),FUN=sum)
allplotsperiode = expand.grid(period=unique(energy$period),plot=unique(energy$plot))
allplotsperiode = merge(allplotsperiode,treatment)

energy = merge(energy,treatment)
energey = merge(energy,allplotsperiode,all=T)
energey$x[is.na(energey$x)]=0

energysummary <- summarySE(energey,measurevar='x',groupvars = c('treatment','period'))

ggplot(energysummary[energysummary$treatment %in% c('CC','EE','EC'),], aes(x=period, y=x, colour=treatment)) + 
  geom_errorbar(aes(ymin=x-se, ymax=x+se), width=.1) +
  geom_line(size=1.5) +
  geom_point() +
  scale_y_continuous(name='Total sm gran energy') +
  scale_x_continuous(name='') +
  geom_vline(xintercept = as.numeric(436),size=1.5) +
  scale_colour_manual(name="Experimental\nTreatment",
                      breaks=c("CC", "EE", "EC"),
                      labels=c("Long-term\nControl\n", "Long-term\nExclosure\n", "Exclosure to\nControl"),
                      values=cbPalette[c(6,4,8)])

# ============================================================================================
# sp richness

