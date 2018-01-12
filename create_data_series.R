library(dplyr)
library(mgcv)
source('data_functions.R')

# Here is some code to use the data_functions.R file to create timeseries of various rodent metrics.

# this function takes a while -- will have to eventually figure out a way to speed it up
data = get_data()




# dipo abundance by plot
dipoN = make_N_data(species = 'Dipos', data)

# dipo abundance; average by treatment
dipoN_avg = avg_by_trt(dipoN)
#write.csv(dipoN_avg,'Dipo_abundance_by_treatment.csv',row.names = F)

# small granivore abundance
sm_gran = make_N_data(species='SmGran', data)

# small granivore abundance; average by treatment
sm_gran_avg = avg_by_trt(sm_gran)

# species richness
sprich = make_speciesrich_data(data)

# species richness; average by treatment
sprich_avg = avg_by_trt(sprich)



## total abundance
#raw_abund = make_data(species='All', start_period=415) %>% trt_data()

## just small heteromyids
#sm_hets = make_data(species='SmH', start_period=415) %>% trt_data()

## just small murids
#sm_mur = make_data(species='SmM', start_period=415) %>% trt_data()

## total rodent mass per plot per period 
#mass_dat = get_mass(start_period = 415,avg=F,metE=F) %>% trt_data()

## average rodent mass per plot per period
#avg_mass = get_mass(start_period = 415,avg=T,metE=F) %>% trt_data()

## total rodent metabolic energy per plot/period
#metE_tot = get_mass(start_period = 415,avg=F,metE=T) %>% trt_data()

## species evenness
#even_dat = get_evenness(start_period = 415) %>% trt_data()



