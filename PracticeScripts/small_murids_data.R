library(dplyr)
library(mgcv)
library(ggplot2)
source('gam_functions.R')
source('data_functions.R')



##### Data for all small granivore analyses
murids = make_data('SmM') %>% trt_data()

CC_mur = murids[[1]]
EC_mur = murids[[2]]
XC_mur = murids[[3]]
CE_mur = murids[[4]]
EE_mur = murids[[5]]

### Plot the data
plot(CC_mur$date,CC_mur$n,xlab='',ylab='small murid abund',main='controls',pch=16,col='red')
avgcc = aggregate(CC_mur$n,by=list(date=CC_mur$date),FUN=mean)
lines(avgcc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EC_mur$date,EC_mur$n,xlab='',ylab='small murid abund',main='krat excl. -> control',pch=16,col='blue')
avgec = aggregate(EC_mur$n,by=list(date=EC_mur$date),FUN=mean)
lines(avgec,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(XC_mur$date,XC_mur$n,xlab='',ylab='small murid abund',main='total excl. -> control',pch=16,col='forestgreen')
avgxc = aggregate(XC_mur$n,by=list(date=XC_mur$date),FUN=mean)
lines(avgxc,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(CE_mur$date,CE_mur$n,xlab='',ylab='small murid abund',main='control -> krat excl.',pch=16,col='purple')
avgce = aggregate(CE_mur$n,by=list(date=CE_mur$date),FUN=mean)
lines(avgce,lwd=2)
abline(v=as.Date('2015-03-15'))

plot(EE_mur$date,EE_mur$n,xlab='',ylab='small murid abund',main='krat excl -> krat excl',pch=16,col='pink')
avgee = aggregate(EE_mur$n,by=list(date=EE_mur$date),FUN=mean)
lines(avgee,lwd=2)
abline(v=as.Date('2015-03-15'))

