#' THis script is to explore whether krats use space differently on control plots vs former krat removals
#' EMC 6/24/19


library(dplyr)

source("FinalAnalysis/population_model_functions.r")

# DATA FILES (PortalData folder was downloaded as part of create_data_series.R script)

# rodent file
rdat <- read.csv('PortalData/Rodents/Portal_rodent.csv', header = TRUE, na.strings = c(""), stringsAsFactors = FALSE)


##########################################################
# DATA PREP
##########################################################
# restrict to only the plots relevant to this project (controls after the switch in 2015)
plotswitchplots = c(4,11,14,17,6,13,18,5,7,24)

rdat_filtered = dplyr::filter(rdat, period>=437, plot %in% plotswitchplots)
spdat = individual_tag_cleanup(sp, rdat_filtered)


# find out how many plots/stakes associated with each individual
dplyr::filter(spdat, group==2)
