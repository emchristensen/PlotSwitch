library(dplyr)
library(mgcv)
library(ggplot2)
source('gam_functions.R')
source('data_functions.R')

# get species richness data
sprich = rodent_abundance(species='All',start_period=400,incomplete=F) %>% species_rich() %>% trt_data()
CC_SR = sprich[[1]]
EC_SR = sprich[[2]]
XC_SR = sprich[[3]]

### Plot the data
plot(CC_SR$date,CC_SR$nsp,xlab='',ylab='# species',main='controls',pch=16,col='red')
avgcc = aggregate(CC_SR$nsp,by=list(date=CC_SR$date),FUN=mean)
lines(avgcc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EC_SR$date,EC_SR$nsp,xlab='',ylab='# species',main='krat excl. -> control',pch=16,col='blue')
avgec = aggregate(EC_SR$nsp,by=list(date=EC_SR$date),FUN=mean)
lines(avgec,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(XC_SR$date,XC_SR$nsp,xlab='',ylab='# species',main='total excl. -> control',pch=16,col='forestgreen')
avgxc = aggregate(XC_SR$nsp,by=list(date=XC_SR$date),FUN=mean)
lines(avgxc,lwd=2)
abline(v=as.Date('2015-03-15'))

##### CC Plots 
m_CC <- gam(nsp ~ s(month, bs = "cc", k = 12) + s(Time), data = CC_SR, family= poisson)
gam_diagnostics(m_CC, "CC no AR")

# create trend info 
CC_trend = make_prediction_gam(CC_SR,m_CC,colname='nsp')

# plot gam results
plot_singleGAM(CC_trend,treatment='CC')
