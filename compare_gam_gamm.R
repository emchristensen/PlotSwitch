# this script is exploratory code, looking at the difference between the gam and gamm functions from mgcv


library(dplyr)
library(RCurl)
library(mgcv)
library(ggplot2)
source('gamm_functions.R')
source('gam_functions.R')
source('data_functions.R')

##### Data for all Dipo analyses
dipo_data = make_dipo_data()
filtered_data = trt_data(dipo_data)
CC = filtered_data[[1]]
EC = filtered_data[[2]]
XC = filtered_data[[3]]

GAM_type = 'uncorrelated errors/gaussian'
Ylab = 'Dipodomys abundance/plot'

fam = Gamma(link=identity)

##### CC Plots -- GAMM
m_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k=20), random=list(plot=~1), data = CC, family=fam, link = identity)
gamm_diagnostics(m_CC, "CC no AR")

# create trend info 
CC_trend = make_prediction_gamm(CC,m_CC)
op <- par(mar = c(5,4,2,2) + 0.1)

# plot gam results
plot_singleGAM(CC_trend, GAM_type, Ylab, "CC")

# with AR1
m1_CC <- gamm(DipoN ~ s(month, bs = "cc", k = 12) + s(Time, k=20),
              data = CC, correlation = corARMA(form = ~ 1|year, p = 1), family = fam)
gamm_diagnostics(m1_CC,"CC AR1")

###### CC Plots -- GAM
m3 = gam(DipoN ~ s(month, bs = "cc", k = 12) + s(Time), data = CC, family=fam, link=identity)
CC_gamtrend = make_prediction_gam(CC,m3,colname='DipoN')
plot_singleGAM(CC_gamtrend, GAM_type, Ylab, "CC")
gam_diagnostics(m3,'')
gam.check(m3)
