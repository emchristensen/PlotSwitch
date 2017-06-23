library(dplyr)
library(mgcv)
library(ggplot2)
source('gam_functions.R')
source('data_functions.R')

# get species richness data
rdat = rodent_abundance(species='All',start_period=400,incomplete=F)

sprich = species_rich(rdat) %>% trt_data()
CC_SR = sprich[[1]]
EC_SR = sprich[[2]]
XC_SR = sprich[[3]]

##### CC Plots 
m_CC <- gam(nsp ~ s(month, bs = "cc", k = 12) + s(Time), data = CC_SR, family= poisson)
gam_diagnostics(m_CC, "CC no AR")

# create trend info 
CC_trend = make_prediction_gam(CC_SR,m_CC,colname='nsp')

# plot gam results
plot_singleGAM(CC_trend,treatment='CC')
