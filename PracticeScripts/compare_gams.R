library(dplyr)
library(RCurl)
library(mgcv)
library(ggplot2)

source('gam_functions.R')
source('data_functions.R')

##### Items for all Dipo analyses
dipo_data = make_dipo_data()
GAM_type = 'uncorrelated errors/gaussian'
Ylab = 'Dipodomys abundance/plot'


##### CC Plots 
CC = dipo_data %>% filter(treatment == "CC") %>% arrange(date)
CC$plot = as.factor(CC$plot)
m_CC_plot <- gam(DipoN ~ plot + s(month, bs = "cc", k = 12) + s(Time), data = CC)
m_CC <- gam(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = CC)

gam_diagnostics(m_CC_plot, "Plot in model")


CC_trend = make_prediction_gam(CC,m_CC,colname='DipoN')
op <- par(mar = c(5,4,2,2) + 0.1)

# plot gam results

plot_singleGAM(CC_trend, GAM_type, Ylab, "CC")
