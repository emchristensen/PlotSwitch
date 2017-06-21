
library(dplyr)

#' Function to tally captures by species, period, and plot
#' This function was made specifically for working with the 2015 plot switch
#' 
#' @param species species desired; either "All" or "Granivore"
#' @param start_period first period number of data desired; default is 130 (1989)
#' 
#' @usage rdat = read.csv('C:/Users/EC/Desktop/git/PortalData/Rodents/Portal_rodent.csv',as.is=T,na.strings = '')
#'        rodents = rodent_abundance(rdat)

rodent_abundance = function(species='All',start_period=130) {
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"
  rdat = read.csv(text=RCurl::getURL(http),as.is=T,na.strings = '')

  if (species=='All') {targetsp = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','RO','SF','SH')}
  if (species=='Granivore') {targetsp = c('BA','DM','DO','DS','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')}
  # filter data by desired species; remove early trapping periods
  rdat_filtered = filter(rdat, species %in% targetsp, period >= start_period)
  
  # aggregate by species, plot, period
  byspecies = aggregate(rdat_filtered$species,by=list(period = rdat_filtered$period, plot = rdat_filtered$plot,species = rdat_filtered$species),FUN=length)
  
  plottype = data.frame(plot = seq(1,24),type_before=c('C','C','E','C','X','E','X','C','C','X','C','C','E','C','E','X','C','E','E','E','E','C','X','X'),
                        type_after=c('X','E','E','C','C','C','C','E','X','X','C','X','C','C','E','X','C','C','E','E','E','E','X','C'),stringsAsFactors = F)
  
  byspecies_merge = merge(byspecies,plottype)
  
  return(byspecies_merge)
}


