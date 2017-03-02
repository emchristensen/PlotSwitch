# LDA analysis (with figures)
# using treatment-level data

# uses data from controls and k-rat exclosures 1989-2015

# LDA =================================
library(topicmodels)
library(ggplot2)
library(lubridate)
library(dplyr)
library(reshape2)
set.seed(20)


# read in latest Rodent data, make a table of species/time
rdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent.csv',stringsAsFactors = F,na.strings = '')
rodents = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PL','PI','PH','PP','RM','RO','RF','SF','SH','SO')

# read in trapping data and find censuses where all control plots were trapped
trapdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent_trapping.csv')
allplots = aggregate(trapdat$Sampled,by=list(period=trapdat$Period),FUN=sum) %>% filter(x>22)

# create table of species abundances: 1989 - 2015 so we can use more k-rat removal plots
absabund = rdat %>% filter(period>380, period< 437, period %in% allplots$period, species %in% rodents, plot %in% c(2,4,8,11,12,14,17,22,3,6,13,15,18,19,20,21)) %>% 
  select(period,plot,species) 
counts = aggregate(absabund$species,by=list(period=absabund$period,plot=absabund$plot,species=absabund$species),FUN=length)
counttable = data.frame(reshape(counts,timevar='species',idvar=c('period','plot'),direction='wide')) %>% arrange(period,plot)
counttable[is.na(counttable)] = 0

# dates of trapping periods
#perdat= read.csv('period_dates_single.csv')
#perdat$date = as.Date(perdat$date,format='%m/%d/%Y')

# treatment of each plot (1989-2015)
treatment = data.frame(plot=seq(1:24),treat = c('C','C','E','C','X','E',
                                                'X','C','C','X','C','C',
                                                'E','C','E','X','C','E',
                                                'E','E','E','C','X','X'))

# sum captures of each species by period and treatment (there are 8 plots of each treatment type)
df = merge(counttable,treatment)
df_melt <- melt(df, id = c("period", "plot", "treat"))
countstreatment = dcast(df_melt, period + treat ~ variable, sum)

# take off period and treat columns for input to LDA model
dat = countstreatment[,-(1:2)]

# LDA models: groups from 2 to 5
nstart = 20 # For the final analysis, maybe do 1000
ldamodel2 = LDA(dat,2,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel3 = LDA(dat,3,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel4 = LDA(dat,4,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel5 = LDA(dat,5,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")
ldamodel6 = LDA(dat,6,control=list(estimate.alpha=F,alpha=.5, nstart = nstart),method="VEM")


# composition of component communities, to identify most important species
structure(round(exp(ldamodel4@beta), 3), dimnames = list(NULL, ldamodel4@terms))


# ===========================================
# figures

#dates = as.Date(perdat$date[1:length(dat[,1])])

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_component_communities = function(ldamodel,k,index,treatment) {
  # function to take output of LDA function and plot time series of component communities
  z = posterior(ldamodel)
  topicresults = cbind(data.frame(z$topics),index)
  
  ldaplot = melt(topicresults,id=c('period','treat'))
  
  ggplot(filter(ldaplot,treat==treatment), aes(x=period,y=value,colour=as.factor(variable))) + 
    geom_point() +
    geom_line(size=1) 
    #geom_point(data=filter(ldaplot,treat=='E')) +
    #geom_line(data=filter(ldaplot,treat=='E'))
    #scale_x_continuous(name='') +
    #scale_y_continuous(name='Percent Similarity') 
    #theme(axis.text=element_text(size=18),
    #      axis.title=element_text(size=18),
    #      legend.text=element_text(size=18),
    #      legend.title=element_text(size=18)) +
    #scale_colour_manual(name="",
    #                    breaks=as.character(seq(k)),
    #                    values=cbPalette[1:k])
  
}

index = countstreatment[,1:2]

# plots
plot_component_communities(ldamodel2,2,index,'E')
plot_component_communities(ldamodel3,3,index,'C')
plot_component_communities(ldamodel4,4,index,'E')
plot_component_communities(ldamodel5,5,index,'E')
plot_component_communities(ldamodel6,6,index,'E')

