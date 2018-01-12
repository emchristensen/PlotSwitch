# GAM practice 5/15/17

# https://www.r-bloggers.com/modelling-seasonal-data-with-gams/


library(dplyr)
library(RCurl)
library(mgcv)
library(ggplot2)
source('gam_functions.R')
source('gamm_functions.R')
source('data_functions.R')

fam = Gamma
ctrl <- list(niterEM = 0, optimMethod="L-BFGS-B", maxIter = 100, msMaxIter = 100)
knots <- list(nMonth = c(0.5, seq(1, 12, length = 10), 12.5))
# with dipo data
dipo_data = make_dipo_data()
CC = dipo_data %>% filter(treatment == "CC") %>% arrange(date)

m_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k= 20),data = CC, family = fam(link=log), knots = knots)

m1_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
           data = CC, correlation = corARMA(form = ~ 1|year, p = 1), family = fam(link=log), knots = knots)
m2_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|year, p = 2), family = fam(link=log), knots = knots)
m3_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|year, p = 3), family = fam(link=log), knots = knots)
m4_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|year, p = 4), family = fam(link=log), knots = knots)
m5_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|year, p = 5),family = fam(link=log), knots = knots)
m6_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = CC, correlation = corARMA(form = ~ 1|year, p = 6),family = fam, knots = knots)

#anova(m_CC$lme, m1_CC$lme)   # This command doesn't work with gamms
AIC(m_CC$lme,m1_CC$lme,m2_CC$lme,m3_CC$lme,m4_CC$lme,m5_CC$lme)
layout(matrix(1:2, ncol = 2))
plot(m2_CC$gam, scale = 0)
acf(resid(m2_CC$lme,type='normalized'), lag.max = 36, main = "ACF")
pacf(resid(m2_CC$lme,type='normalized'), lag.max = 36, main = "pACF")
layout(1)

CC2_trend = make_prediction_gamm(CC,m2_CC)
plot(DipoN  ~ date, data = CC, type = "p", 
     ylab = '')
lines(p_raw ~ date, data = CC2_trend, col = "blue")

# EC plots
EC = dipo_data %>% filter(treatment == "EC") %>% arrange(date)
m_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time),data = EC)
m1_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|year, p = 1))
m2_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|year, p = 2))
m3_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|year, p = 3))
m4_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|year, p = 4))
m5_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|year, p = 5),
              control = ctrl)
m7_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC, correlation = corARMA(form = ~ 1|year, p = 7))

anova(m_EC$lme, m1_EC$lme, m2_EC$lme, m3_EC$lme, m4_EC$lme,m7_EC$lme)
gamm_diagnostics(m1_CC,'')

# ===========================================================================================
# EC plots starting at flip
EC2 = read.csv('EC_startatswitch.csv')  # this file contains interpolations to estimate missing data from skipped months
EC2$date = as.Date(EC2$date,format='%m/%d/%Y')
m_EC2 <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time),data = EC2)
m1_EC2 <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC2, correlation = corARMA(form = ~ 1|year, p = 1))
m2_EC2 <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC2, correlation = corARMA(form = ~ 1|year, p = 2))
m3_EC2 <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC2, correlation = corARMA(form = ~ 1|year, p = 3))
m4_EC2 <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC2, correlation = corARMA(form = ~ 1|year, p = 4))
m5_EC2 <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC2, correlation = corARMA(form = ~ 1|year, p = 5))
m6_EC2 <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k = 20),
              data = EC2, correlation = corARMA(form = ~ 1|year, p = 6))

anova(m_EC2$lme, m1_EC2$lme, m2_EC2$lme, m3_EC2$lme, m4_EC2$lme,m5_EC2$lme)
gamm_diagnostics(m_EC2,'')

EC2_trend = make_prediction_gamm(EC2,m_EC2)
plot_singleGAM(EC2_trend, 'gaussian', 'n dipo', "EC")

XC2 = read.csv('XC_startatswitch.csv')  # this file contains interpolations to estimate missing data from skipped months
XC2$date = as.Date(XC2$date,format='%m/%d/%Y')
m_XC2 <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time),data = XC2)
gamm_diagnostics(m_XC2,'')

XC2_trend = make_prediction_gamm(XC2,m_XC2)
plot_singleGAM(XC2_trend, 'gaussian', 'n dipo', "XC")

# =============================================================
# http://multithreaded.stitchfix.com/blog/2015/07/30/gam/
# Kim Larsen

set.seed(3)
x <- seq(0,2*pi,0.1)
z <- sin(x)
y <- z + rnorm(mean=0, sd=0.5*sd(z), n=length(x))
d <- cbind.data.frame(x,y,z)
