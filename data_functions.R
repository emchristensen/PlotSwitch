library(dplyr)


#' Function to tally captures by species, period, and plot
#' This function was made specifically for working with the 2015 plot switch
#' 
#' @param species species desired; either "All" or "Granivore"
#' @param start_period first period number of data desired; default is 130 (1989)
#' @param incomplete T/F, wheter or not to include incomplete censuses
#' 
#' @return data frame with columns 'plot', 'period', 'species', 'x' (abundance), 'date','month','Year','Time'

rodent_abundance = function(species='All',start_period=130,incomplete=F) {
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"
  rdat = read.csv(text=RCurl::getURL(http),as.is=T,na.strings = '')

  if (species=='All') {targetsp = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','RO','SF','SH')}
  if (species=='Granivore') {targetsp = c('BA','DM','DO','DS','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')}
  
  # Get trapping history in order to find and remove incomplete censuses
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"
  trapping = read.csv(text=RCurl::getURL(http)) 
  trapping$date = as.Date(paste(trapping$Year,trapping$Month,trapping$Day,sep='-'))
  plotstrapped = aggregate(trapping$Sampled,by=list(period=trapping$Period),FUN=sum)
  fullcensus = plotstrapped[plotstrapped$x>=21,] #I'm using 21 instead of 24 because  I know period 457 only trapped 21 plots, but the skipped plots aren't used in this project
  
  # filter data by desired species; remove early trapping periods
  rdat_filtered = dplyr::filter(rdat, species %in% targetsp, period >= start_period)
  
  # if desired, remove incompete trapping periods
  if (incomplete==F) {
    rdat_filtered = filter(rdat_filtered, period %in% fullcensus$period)
  }
  
  # aggregate by species, plot, period
  byspecies = aggregate(rdat_filtered$species,by=list(period = rdat_filtered$period, plot = rdat_filtered$plot,species = rdat_filtered$species),FUN=length)
  
  # adding sampling date to the dipo table
  #first_date = trapping %>% select(Period,date) %>% group_by(Period) %>% summarise_each(funs(min), date)
  first_date = aggregate(trapping$date,by=list(Period=trapping$Period),FUN=min)
  names(first_date) = c('Period','date')
  byspecies_date = left_join(byspecies, first_date, by=c('period'='Period'))
  byspecies_date = byspecies_date %>% mutate(month = as.numeric(format(date, "%m")),
                                 Year = as.numeric(format(date, "%Y")),
                                 Time = as.numeric(date) / 1000)
  
  return(byspecies_date)
}


#' @title Make Dipo Data
#' 
#' @description make a data frame for running GAM models
#' 
#' @param
#' 
#' @return data frame with columns:
#'            plot
#'            period
#'            treatment (two letters representing before/after switch: C = control, E = krat exclosure, X = total rodent removal)
#'            DipoN (number of DM, DO, DS summed)
#'            date 
#'            month
#'            Year
#'            Time

make_dipo_data = function(){
  # Create species abundances by plot by period
  byspecies = rodent_abundance(species='Granivore',incomplete=F)
  
  # number of dipos per plot
  d = byspecies %>% 
    filter(period>414, species %in% c('DM','DO','DS'))
  dipos = aggregate(d$x, by=list(period=d$period, plot=d$plot, date=d$date, month=d$month, Year=d$Year, Time=d$Time), FUN=sum)
  # data frame of all plots in all periods
  allplotsperiod = expand.grid(period=unique(dipos$period), plot=unique(dipos$plot))
  #attach date columns according to period
  dateinfo = unique(dipos[,c('period','date','month','Year','Time')])
  allplotsperiod = merge(allplotsperiod,dateinfo)
  # attach treatment column according to plot number
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  allplotsperiod = merge(allplotsperiod,treatment)
  # merge capture data with data frame of all plots and all periods, and fill in empty data with zeros
  dipos = merge(allplotsperiod,dipos,all=T)
  dipos$x[is.na(dipos$x)] = 0
  # change column name
  dipos_gam =  rename(dipos,DipoN=x)
  # put data in chronological order
  dipos_gam = dipos_gam[order(dipos_gam$period),]
  
  return(dipos_gam)
}


#' @title Treatment Data
#' 
#' @description takes tallies of dipos from make_dipo_data and separates by treatment
#' 
#' @param data dataframe output from make_dipo_data
#' 
#' @return a list containing three data frames, containing dipo data from 3 treatment types

trt_data = function(data){
  data$plot = as.factor(data$plot)
  CC = data %>% filter(treatment == "CC") %>% arrange(date)
  EC = data %>% filter(treatment == "EC") %>% arrange(date)
  XC = data %>% filter(treatment == "XC") %>% arrange(date)
  return(list(CC,EC,XC))
}



#' @title Species richness
#' 
#' @description get species richness at period and plot level
#' 
#' @param rdat rodent data; column period, plot, species
#' 
#' @return data frame with period, plot, nsp (number of species)

species_rich = function(rdat) {
  
  # count number of species in rodent data
  richness = aggregate(rdat$species,
                       by=list(period=rdat$period,plot=rdat$plot,date=rdat$date,month=rdat$month,Year=rdat$Year,Time=rdat$Time),
                       FUN=unique)
  for (n in 1:length(richness$period)) {
    richness$nsp[n] = length(unlist(richness$x[n]))
  }
  
  # make sure there are zeros where plot was trapped but nothing caught; add column for treatment type
  allplotsperiod = expand.grid(period=unique(sprich$period), plot=unique(sprich$plot))
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  allplotsperiod = merge(allplotsperiod,treatment)
  sprich = merge(allplotsperiod,richness,all=T)
  sprich$nsp[is.na(sprich$nsp)] = 0
  
  return(dplyr::select(sprich,period,plot,nsp,treatment,date,month,Year,Time))
}
