library(dplyr)
library(mgcv)
library(ggplot2)
source('gam_functions.R')
source('data_functions.R')

# Here is some code to use the data_functions.R file to create timeseries of various rodent metrics.
# In all cases, the column of interest is called 'n', which will hopefully make it easier to plug
# new data into existing code.  


# dipo abundance
dipo_dat = make_data(species='Dipos',start_period=415) %>% trt_data()
# or
dipo_dat = make_dipo_data()

# total abundance
raw_abund = make_data(species='All', start_period=415) %>% trt_data()

# get small granivore abundance
sm_gran = make_data(species='SmGran', start_period=415) %>% trt_data()

# just small heteromyids
sm_hets = make_data(species='SmH', start_period=415) %>% trt_data()

# just small murids
sm_mur = make_data(species='SmM', start_period=415) %>% trt_data()

# species richness
sprich = rodent_abundance(start_period=415,incomplete=F) %>% species_rich() %>% trt_data()

# total rodent mass per plot per period 
mass_dat = get_mass(start_period = 415,avg=F,metE=F) %>% trt_data()

# average rodent mass per plot per period
avg_mass = get_mass(start_period = 415,avg=T,metE=F) %>% trt_data()

# total rodent metabolic energy per plot/period
metE_tot = get_mass(start_period = 415,avg=F,metE=T) %>% trt_data()

