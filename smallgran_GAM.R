library(dplyr)
library(mgcv)
library(ggplot2)
source('gam_functions.R')
source('data_functions.R')



##### Data for all small granivore analyses
smgran = make_data('SmGran') %>% trt_data()

CC_small = smgran[[1]]
EC_small = smgran[[2]]
XC_small = smgran[[3]]
CE_small = smgran[[4]]
EE_small = smgran[[5]]


### Plot the data
plot(CC_small$date,CC_small$n,xlab='',ylab='small gran abund',main='controls',pch=16,col='red')
avgcc = aggregate(CC_small$n,by=list(date=CC_small$date),FUN=mean)
lines(avgcc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EC_small$date,EC_small$n,xlab='',ylab='small gran abund',main='krat excl. -> control',pch=16,col='blue')
avgec = aggregate(EC_small$n,by=list(date=EC_small$date),FUN=mean)
lines(avgec,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(XC_small$date,XC_small$n,xlab='',ylab='small gran abund',main='total excl. -> control',pch=16,col='forestgreen')
avgxc = aggregate(XC_small$n,by=list(date=XC_small$date),FUN=mean)
lines(avgxc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(CE_small$date,CE_small$n,xlab='',ylab='small gran abund',main='control -> krat excl.',pch=16,col='purple')
avgce = aggregate(CE_small$n,by=list(date=CE_small$date),FUN=mean)
lines(avgce,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EE_small$date,EE_small$n,xlab='',ylab='small gran abund',main='krat excl -> krat excl',pch=16,col='pink')
avgee = aggregate(EE_small$n,by=list(date=EE_small$date),FUN=mean)
lines(avgee,lwd=2)
abline(v=as.Date('2015-03-15'))
