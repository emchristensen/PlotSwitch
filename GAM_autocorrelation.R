# GAM practice 5/15/17

# https://www.r-bloggers.com/modelling-seasonal-data-with-gams/


library(dplyr)
library(RCurl)
library(mgcv)
library(ggplot2)
source('gam_functions.R')


# with dipo data
dipo_data = make_dipo_data()
CC = dipo_data %>% filter(treatment == "CC") %>% arrange(date)
m_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time),data = CC)
m1_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
           data = CC, correlation = corARMA(form = ~ 1|Year, p = 1))
m2_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|Year, p = 2))
m3_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|Year, p = 3))
m4_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|Year, p = 4))
m5_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|Year, p = 5))
m6_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|Year, p = 6))

anova(m_CC$lme, m1_CC$lme, m2_CC$lme, m3_CC$lme, m4_CC$lme, m5_CC$lme, m6_CC$lme)

layout(matrix(1:2, ncol = 2))
plot(m5_CC$gam, scale = 0)
acf(resid(m5_CC$lme,type='normalized'), lag.max = 36, main = "ACF")
pacf(resid(m5_CC$lme,type='normalized'), lag.max = 36, main = "pACF")
layout(1)

CC5_trend = make_prediction(CC,m5_CC)
plot(DipoN  ~ date, data = CC, type = "p", 
     ylab = ylab)
lines(p_raw ~ date, data = CC5_trend, col = "blue")

# EC plots
EC = dipo_data %>% filter(treatment == "EC") %>% arrange(date)
m_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time),data = EC)
m1_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|Year, p = 1),
              control = ctrl)
m2_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|Year, p = 2),
              control = ctrl)
m3_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|Year, p = 3),
              control = ctrl)
m4_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|Year, p = 4),
              control = ctrl)
m5_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|Year, p = 5),
              control = ctrl)
m7_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|Year, p = 7),
              control = ctrl)

anova(m_EC$lme, m1_EC$lme, m2_EC$lme, m3_EC$lme, m4_EC$lme, m5_EC$lme, m6_EC$lme)
gam_diagnostics(m7_CC,'EC AR5')


# =============================================================
# http://multithreaded.stitchfix.com/blog/2015/07/30/gam/
# Kim Larsen

set.seed(3)
x <- seq(0,2*pi,0.1)
z <- sin(x)
y <- z + rnorm(mean=0, sd=0.5*sd(z), n=length(x))
d <- cbind.data.frame(x,y,z)
