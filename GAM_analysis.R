library(dplyr)
library(RCurl)
library(mgcv)

source('gam_functions.R')

ylab = "# of Dipos"

dipo_data = make_dipo_data()
  
# Seasonal GAM on CC plots only, no autocorrelation correction
CC = dipo_data %>% filter(treatment == "CC") %>% arrange(date)
m_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = CC)

gam_diagnostics(m_CC, "CC no AR")

# plot trend on data
CC_trend = make_prediction(CC,m_CC)
op <- par(mar = c(5,4,2,2) + 0.1)

plot(DipoN  ~ date, data = CC, type = "p", 
     ylab = ylab)
abline(v=as.Date("2015-03-15", format="%Y-%m-%d"))
lines(p_raw  ~ date, data = CC_trend, col = "black")
legend("topleft", legend = c("Uncorrelated Errors"), lty = 1, lwd = c(1,1,1))
par(op)


######  EC PLOTS ###############
EC = dipo_data %>% filter(treatment == "EC") %>% arrange(date)
plot(DipoN ~ date, data = EC, type = "p", ylab = ylab)
abline(v=16.52)

# Seasonal GAM on EC plots only, no autocorrelation correction
m_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = EC)
gam_diagnostics(m_EC, "EC no AR")

EC_trend = make_prediction(EC,m_EC)

op <- par(mar = c(5,4,2,2) + 0.1)

plot(DipoN  ~ date, data = EC, type = "p", 
     ylab = ylab)
abline(v=as.Date("2015-03-15", format="%Y-%m-%d"))
lines(p_raw  ~ date, data = EC_trend, col = "black")
lines(p_raw ~ date, data = CC_trend, col = "blue")
legend("topleft", legend = c("Uncorrelated Errors"), lty = 1, lwd = c(1,1,1))
par(op)

######  XC PLOTS ###############
XC = dipo_data %>% filter(treatment == "XC") %>% arrange(date)
plot(DipoN ~ date, data = XC, type = "p", ylab = ylab)
abline(v=16.52)

# Seasonal GAM on EC plots only, no autocorrelation correction
m_XC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = XC)
gam_diagnostics(m_XC, "XC no AR")

XC_trend = make_prediction(XC,m_XC)

op <- par(mar = c(5,4,2,2) + 0.1)

plot(DipoN  ~ date, data = XC, type = "p", 
     ylab = ylab)
abline(v=as.Date("2015-03-15", format="%Y-%m-%d"))
lines(p_raw  ~ date, data = XC_trend, col = "black")
lines(p_raw ~ date, data = CC_trend, col = "blue")
lines(p_raw ~ date, data = EC_trend, col = 'green')
legend("topleft", legend = c("Uncorrelated Errors"), lty = 1, lwd = c(1,1,1))
par(op)
