############ Running GAMs on Plot Switch ########
#
# Current state: technically fits GAMs to each treatment
#                Dipo data only
# TO DO: 1) Add negative binomial model
#        2) Add autocorrelation residual structure
#        3) Plot confidence intervals on trend fits
#################################################

library(dplyr)
library(RCurl)
library(mgcv)
library(ggplot2)

source('gam_functions.R')

ylab = "# of Dipos"

dipo_data = make_dipo_data()
  
##### CC Plots 
CC = dipo_data %>% filter(treatment == "CC") %>% arrange(date)
CC$plot = as.factor(CC$plot)
m_CC_plot <- gamm(DipoN ~ plot + s(month, bs = "cc", k = 12) + s(Time), data = CC)
m_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = CC)

# plot trend on data
CC_trend = make_prediction(CC,m_CC)
op <- par(mar = c(5,4,2,2) + 0.1)

transition = as.Date("2015-03-15", format="%Y-%m-%d")
plot_CC = ggplot(aes(x=date, y=p_raw), data = CC_trend) +
  geom_ribbon(aes(ymin=lower, ymax=upper), data= CC_trend,fill='gray90') +
  geom_line(color='blue') +
  geom_vline(xintercept =  as.numeric(transition)) +
  ggtitle("Dipodomys response to plot flip (uncorrelated errors & gaussian)") +
  xlab("Date") + ylab("Dipodomys abundance per plot") +
  theme_classic()

plot_CC
              
######  EC PLOTS ###############
EC = dipo_data %>% filter(treatment == "EC") %>% arrange(date)
plot(DipoN ~ date, data = EC, type = "p", ylab = ylab)
abline(v=16.52)

# Seasonal GAM on EC plots 
m_EC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = EC)
gam_diagnostics(m_EC, "EC no AR")

EC_trend = make_prediction(EC,m_EC)

op <- par(mar = c(5,4,2,2) + 0.1)

######  XC PLOTS ###############
XC = dipo_data %>% filter(treatment == "XC") %>% arrange(date)
plot(DipoN ~ date, data = XC, type = "p", ylab = ylab)
abline(v=16.52)

# Seasonal GAM on XC plots 
m_XC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = XC)
gam_diagnostics(m_XC, "XC no AR")

XC_trend = make_prediction(XC,m_XC)

op <- par(mar = c(5,4,2,2) + 0.1)

######## Plot all plots together, with CI ################
ggplot(aes(x=date, y=p_raw), data = CC_trend) +
  geom_ribbon(aes(ymin=lower, ymax=upper), data= CC_trend,fill='gray90') +
  geom_line(color='blue') +
  geom_ribbon(aes(ymin=lower, ymax=upper), data= EC_trend,fill='gray90') +
  geom_line(color='green', data = EC_trend) +
  geom_ribbon(aes(ymin=lower, ymax=upper), data= XC_trend,fill='gray90') +
  geom_line(color='red', data = XC_trend) +
  geom_vline(xintercept =  as.numeric(transition)) +
  ggtitle("Dipodomys response to plot flip (uncorrelated errors & gaussian)") +
  xlab("Date") + ylab("Dipodomys abundance per plot") +
  theme_classic()


