#' THis script is to explore whether krats use space differently on control plots vs former krat removals
#' EMC 6/24/19

library(ggplot2)
library(dplyr)

source("FinalAnalysis/population_model_functions.r")


# DATA FILES (PortalData folder was downloaded as part of create_data_series.R script)

# rodent file
rdat <- read.csv('PortalData/Rodents/Portal_rodent.csv', header = TRUE, na.strings = c(""), stringsAsFactors = FALSE)

pdat <- data.frame(plot = 1:24, treatment = c('CC','CE','EE','CC','XC','EC',
                                              'XC','CE','CX','XX','CC','CX',
                                              'EC','CC','EE','XX','CC','EC',
                                              'EE','EE','EE','CE','XX','XC'))
##########################################################
# DATA PREP
##########################################################
# restrict to only the plots relevant to this project (controls after the switch in 2015)
plotswitchplots = c(4,11,14,17,6,13,18,5,7,24)

dat_filtered = dplyr::filter(rdat, period>=437, plot %in% plotswitchplots)
sp = 'DM'
spdat = individual_tag_cleanup(sp, rdat_filtered)


# find out how many plots/stakes associated with each individual
dplyr::filter(spdat, group==2)
number_plots_stakes = spdat %>% group_by(group, plot) %>%
  summarise(nstakes = n_distinct(stake),
            ncaptures = n_distinct(recordID)) %>% 
  merge(pdat, by='plot') %>%
  dplyr::filter(ncaptures>1)

# find individuals that are found in multiple plots
plots_per_individual <- spdat %>% group_by(group) %>%
  summarise(nplots = n_distinct(plot))


ggplot(number_plots_stakes, aes(x=ncaptures, y=nstakes, color=treatment)) + geom_jitter()

#########################################################
# raster plot
#########################################################

source('FinalAnalysis/raster_plot_species_by_stake.R')

path = './'
group_or_individual = 'DM'
min_period = 437
raster_of_plot_captures(path, group_or_individual, min_period)