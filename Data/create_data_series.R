# This script uses the data_functions.R file to create timeseries of various rodent metrics.

source('Data/data_functions.R')

## install older version of portalr if necessary
# require(devtools)
#install_version('portalr', version = '0.2.6')

# set location to put PortalData
data_folder = '.'
# Download the PortalData repo, version 1.53.0 (this will take several minutes, but only needs to be done once)
portalr::download_observations(path = data_folder, version = '1.53.0')
# get most recent version of data if wanted
#portalr::download_observations(path = data_folder, version = 'latest')


# Plot treatments we are interested in for this project
treatments = c('CC','EC','XC')

# get abundance data
dat = get_data(startdate = "2013-03-11",
               min_num_plots = 21,
               treatments, 
               data_folder = data_folder)



# dipo abundance by plot
dipoN = make_N_data(species = 'Dipos', dat)
write.csv(dipoN,'Data/Dipo_counts.csv',row.names=F)


# species richness by plot -- not used in final analysis
#sprich = make_speciesrich_data(dat)
#write.csv(sprich,'Data/SpeciesRichness.csv',row.names=F)


# total metabolic energy by plot (all granivores combined)
total_energy = get_community_energy(path = data_folder,
                                    startdate = "2013-03-11",
                                    min_num_plots = 21,
                                    species='Granivore')
total_energy = dplyr::filter(total_energy,treatment %in% c('CC','EC','XC'))
write.csv(total_energy,'Data/TotalCommunityEnergy.csv',row.names=F)


# abundance of small granivores.
smgran = make_N_data(species='SmGran',dat)
write.csv(smgran,'Data/SmallGranivores.csv',row.names=F)




# plants
winterannuals = make_plant_table(path = data_folder,
                                 selected_plots=1:24,
                                 plant_community='Winter Annuals',
                                 summer_winter = 'winter',
                                 threshold = 0.33)
summerannuals = make_plant_table(path = data_folder,
                                 selected_plots=1:24,
                                 plant_community='Summer Annuals',
                                 summer_winter = 'summer',
                                 threshold = 0.33)
treat_table = make_treatment_table()
 
dat.winter = merge(winterannuals,treat_table, by='plot') %>%
  dplyr::filter(year %in% c(2008,2012,2013,2014,2015))
dat.summer = merge(summerannuals,treat_table, by='plot') %>%
  dplyr::filter(year %in% c(2008,2011,2014))

write.csv(dat.winter,'PlantAnalysis/WinterAnnualTreatments.csv',row.names=F)
write.csv(dat.summer,'PlantAnalysis/SummerAnnualTreatments.csv',row.names=F)


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



