source('data_functions.R')

# Code to use the data_functions.R file to create timeseries of various rodent metrics.

treatments = c('CC','EC','XC')

# this function takes a while -- will have to eventually figure out a way to speed it up
dat = get_data(startdate = "2013-03-11",
               min_num_plots = 21,
               treatments)



# dipo abundance by plot
dipoN = make_N_data(species = 'Dipos', dat)
write.csv(dipoN,'Dipo_counts.csv',row.names=F)


# species richness by plot
sprich = make_speciesrich_data(dat)
write.csv(sprich,'SpeciesRichness.csv',row.names=F)


# total metabolic energy by plot (all granivores combined)
total_energy = get_community_energy(startdate = "2013-03-11",
                                    min_num_plots = 21,
                                    species='Granivore')
total_energy = dplyr::filter(total_energy,treatment %in% c('CC','EC','XC'))
write.csv(total_energy,'TotalCommunityEnergy.csv',row.names=F)


# abund of all small granivores (not just the 5 used in Heske et al 1994) [PP, PE, PF, PM, RM].
smgran = make_N_data(species='SmGran',dat)
write.csv(smgran,'SmallGranivores.csv',row.names=F)









# ==============================================================================
# other, currently unused metrics
#

# total energy of small granivores (all, not just the 5)
#smgran_energy = get_community_energy(startdate = "2013-03-11",min_num_plots = 21,species='SmGran')
#write.csv(smgran_energy,'SmallGranivoreEnergy.csv',row.names=F)

#quarterly = avg_3month(smgran)

# dipo abundance; average by treatment
#dipoN_avg = avg_by_trt(dipoN)
#dipoN_avg = rename(dipoN_avg, dipos=n)
#write.csv(dipoN_avg,'Dipo_abundance_by_treatment.csv',row.names = F)

# DM abundance by plot
#DM = make_N_data(species = 'DM', data)
#DM = rename(DM,dipos=n)
#write.csv(DM,'DM_counts.csv',row.names=F)

# small granivore abundance
#sm_gran = make_N_data(species='SmGran', data)

# small granivore abundance; average by treatment
#sm_gran_avg = avg_by_trt(sm_gran)

# species richness; average by treatment
#sprich_avg = avg_by_trt(sprich)

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



