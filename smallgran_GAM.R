library(dplyr)
library(mgcv)
library(ggplot2)
source('gam_functions.R')
source('data_functions.R')



##### Data for all small granivore analyses
dat = rodent_abundance(species='SmGran',start_period=400,incomplete=F)
totalsg = aggregate(dat$x,by=list(period=dat$period,plot=dat$plot,date=dat$date,month=dat$month,Year=dat$Year,Time=dat$Time),FUN=sum)
# attach treatment column according to plot number
treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                     'XC','EC','XC','CE',
                                     'CX','XX','CC','CX',
                                     'EC','CC','EE','XX',
                                     'CC','EC','EE','EE',
                                     'EE','CE','XX','XC'),plot=seq(1,24))
totalsg= merge(totalsg,treatment)
filtered_data = trt_data(totalsg)
CC_small = filtered_data[[1]]
EC_small = filtered_data[[2]]
XC_small = filtered_data[[3]]
CE_small = filtered_data[[4]]
EE_small = filtered_data[[5]]

### Plot the data
plot(CC_small$date,CC_small$x,xlab='',ylab='small gran abund',main='controls',pch=16,col='red')
avgcc = aggregate(CC_small$x,by=list(date=CC_small$date),FUN=mean)
lines(avgcc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EC_small$date,EC_small$x,xlab='',ylab='small gran abund',main='krat excl. -> control',pch=16,col='blue')
avgec = aggregate(EC_small$x,by=list(date=EC_small$date),FUN=mean)
lines(avgec,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(XC_small$date,XC_small$x,xlab='',ylab='small gran abund',main='total excl. -> control',pch=16,col='forestgreen')
avgxc = aggregate(XC_small$x,by=list(date=XC_small$date),FUN=mean)
lines(avgxc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(CE_small$date,CE_small$x,xlab='',ylab='small gran abund',main='control -> krat excl.',pch=16,col='purple')
avgce = aggregate(CE_small$x,by=list(date=CE_small$date),FUN=mean)
lines(avgce,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EE_small$date,EE_small$x,xlab='',ylab='small gran abund',main='krat excl -> krat excl',pch=16,col='pink')
avgee = aggregate(EE_small$x,by=list(date=EE_small$date),FUN=mean)
lines(avgee,lwd=2)
abline(v=as.Date('2015-03-15'))
