############ Running GAMs on Plot Switch ########
#
# Current state: technically fits GAMs to each treatment
#                Dipo data only
# TO DO: 1) Find acceptable means to evaluate/choose models
#           - autocorrelation?  family?
#           - can fit different models to different treatment types
#
#################################################

library(dplyr)
library(RCurl)
library(mgcv)
library(ggplot2)
source('gam_functions.R')


GAM_type = 'uncorrelated errors/gaussian'
Ylab = 'Dipodomys abundance/plot'

##### Data for all Dipo analyses
dipo_data = make_dipo_data()
filtered_data = trt_data(dipo_data)
CC = filtered_data[[1]]
EC = filtered_data[[2]]
XC = filtered_data[[3]]

##### CC Plots 
m_CC <- gam(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = CC, family=Gamma)
gam_diagnostics(m_CC, "CC no AR")

# create trend info 
CC_trend = make_prediction_gam(CC,m_CC)
op <- par(mar = c(5,4,2,2) + 0.1)

# plot gam results
plot_singleGAM(CC_trend, GAM_type, Ylab, "CC")

######  EC PLOTS ###############
plot(DipoN ~ date, data = EC, type = "p", ylab = ylab)


# Seasonal GAM on EC plots 
m_EC <- gam(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = EC)
gam_diagnostics(m_EC, "EC no AR")

EC_trend = make_prediction_gam(EC,m_EC)

op <- par(mar = c(5,4,2,2) + 0.1)

plot_singleGAM(EC_trend, GAM_type, Ylab, "EC")
abline(v=16.52)
######  XC PLOTS ###############

plot(DipoN ~ date, data = XC, type = "p", ylab = ylab)

# Seasonal GAM on XC plots 
m_XC <- gam(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = XC, family=poisson)
gam_diagnostics(m_XC, "XC no AR")

XC_trend = make_prediction_gam(XC,m_XC)
op <- par(mar = c(5,4,2,2) + 0.1)

plot_singleGAM(XC_trend, GAM_type, Ylab, "XC")
abline(v=16.52)
######## Plot all plots together, with CI ################
transition = as.Date("2015-03-15", format="%Y-%m-%d")
ggplot(aes(x=date, y=p_raw), data = CC_trend) +
  geom_ribbon(aes(ymin=lower, ymax=upper), data= CC_trend,fill='gray90',alpha=.5) +
  geom_line(color='red') +
  geom_ribbon(aes(ymin=lower, ymax=upper), data= EC_trend,fill='gray90',alpha=.5) +
  geom_line(color='green', data = EC_trend) +
  geom_ribbon(aes(ymin=lower, ymax=upper), data= XC_trend,fill='gray90',alpha=.5) +
  geom_line(color='blue', data = XC_trend) +
  geom_vline(xintercept =  as.numeric(as.Date('2015-03-20'))) +
  ggtitle("Dipodomys response to plot flip (uncorrelated errors & gaussian)") +
  xlab("Date") + ylab("Dipodomys abundance per plot") +
  theme_classic() +
  geom_point(aes(x=date,y=DipoN,color=treatment),data=dipo_data)

####################################
#GAM evaluation
gam.check(m_CC) # gives qqplots, residulas
plot(residuals(m_CC)) # plot residuals over time
plot(m_EC,residuals=T,pch=19)
concurvity(m_CC)
