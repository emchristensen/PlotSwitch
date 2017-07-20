library(dplyr)
library(mgcv)
library(ggplot2)
source('gam_functions.R')
source('data_functions.R')



##### Data for small heteromyid analyses
hets = make_data('SmH') %>% trt_data()

CC_het = hets[[1]]
EC_het = hets[[2]]
XC_het = hets[[3]]
CE_het = hets[[4]]
EE_het = hets[[5]]

### Plot the data
plot(CC_het$date,CC_het$n,xlab='',ylab='small heteromyid abund',main='controls',pch=16,col='red')
avgcc = aggregate(CC_het$n,by=list(date=CC_het$date),FUN=mean)
lines(avgcc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EC_het$date,EC_het$n,xlab='',ylab='small heteromyid abund',main='krat excl. -> control',pch=16,col='blue')
avgec = aggregate(EC_het$n,by=list(date=EC_het$date),FUN=mean)
lines(avgec,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(XC_het$date,XC_het$n,xlab='',ylab='small heteromyid abund',main='total excl. -> control',pch=16,col='forestgreen')
avgxc = aggregate(XC_het$n,by=list(date=XC_het$date),FUN=mean)
lines(avgxc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(CE_het$date,CE_het$n,xlab='',ylab='small heteromyid abund',main='control -> krat excl.',pch=16,col='purple')
avgce = aggregate(CE_het$n,by=list(date=CE_het$date),FUN=mean)
lines(avgce,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EE_het$date,EE_het$n,xlab='',ylab='small heteromyid abund',main='krat excl -> krat excl',pch=16,col='pink')
avgee = aggregate(EE_het$n,by=list(date=EE_het$date),FUN=mean)
lines(avgee,lwd=2)
abline(v=as.Date('2015-03-15'))
