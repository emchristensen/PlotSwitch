
library(dplyr)

# ======================================================================================
# Function to tally captures by species, period, and plot
# Parameter: species='All' for all rodent species
#                    'Granivore' for just granivores

# This function is specifically for working with the 2015 plot switch; it excludes periods
# before 1988 (periods 1-129)

rodent_abundance = function(species='All') {
  library(RCurl)
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"
  rdat = read.csv(text=getURL(http),as.is=T,na.strings = '')
  
  # filter data by desired species; remove early trapping periods; exclude "Removed" animals (R in note5)
  if (species=='All') {targetsp = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','RO','SF','SH')}
  if (species=='Granivore') {targetsp = c('BA','DM','DO','DS','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')}
  
  rdat = filter(rdat, species %in% targetsp, period > 129)
  rdat = rdat[(is.na(rdat$note5) | rdat$note5 %in% c('E','D')),]
  
  # aggregate by species, plot, period
  byspecies = aggregate(rdat$species,by=list(period = rdat$period, plot = rdat$plot,species = rdat$species),FUN=length)
  
  plottype = data.frame(plot = seq(1,24),type_before=c('C','C','E','C','X','E','X','C','C','X','C','C','E','C','E','X','C','E','E','E','E','C','X','X'),
                        type_after=c('X','E','E','C','C','C','C','E','X','X','C','X','C','C','E','X','C','C','E','E','E','E','X','C'),stringsAsFactors = F)
  
  byspecies_merge = merge(byspecies,plottype)
  
  return(byspecies_merge)
}

