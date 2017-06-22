library(dplyr)


#' Function to tally captures by species, period, and plot
#' This function was made specifically for working with the 2015 plot switch
#' 
#' @param species species desired; either "All" or "Granivore"
#' @param start_period first period number of data desired; default is 130 (1989)
#' 
#' @return data frame with columns 'plot', 'period', 'species', 'x' (abundance)

rodent_abundance = function(species='All',start_period=130) {
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"
  rdat = read.csv(text=RCurl::getURL(http),as.is=T,na.strings = '')

  if (species=='All') {targetsp = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','RO','SF','SH')}
  if (species=='Granivore') {targetsp = c('BA','DM','DO','DS','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')}
  
  # filter data by desired species; remove early trapping periods
  rdat_filtered = filter(rdat, species %in% targetsp, period >= start_period)
  
  # aggregate by species, plot, period
  byspecies = aggregate(rdat_filtered$species,by=list(period = rdat_filtered$period, plot = rdat_filtered$plot,species = rdat_filtered$species),FUN=length)
  
  return(byspecies)
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
  byspecies = rodent_abundance('Granivore')
  
  # Get trapping history in order to find and remove incomplete censuses
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"
  trapping = read.csv(text=RCurl::getURL(http)) 
  trapping$date = as.Date(paste(trapping$Year,trapping$Month,trapping$Day,sep='-'))
  plotstrapped = aggregate(trapping$Sampled,by=list(period=trapping$Period),FUN=sum)
  fullcensus = plotstrapped[plotstrapped$x>20,]
  
  # number of dipos per plot
  d = byspecies %>% 
    filter(period>414, period %in% fullcensus$period, species %in% c('DM','DO','DS')) %>% 
    select(plot,period,species,x)
  dipos = aggregate(d$x, by=list(period=d$period, plot=d$plot), FUN=sum)
  allplotsperiod = expand.grid(period=unique(dipos$period), plot=unique(dipos$plot))
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  allplotsperiod = merge(allplotsperiod,treatment)
  dipos = merge(allplotsperiod,dipos,all=T)
  dipos$x[is.na(dipos$x)] = 0
  
  # adding sampling date to the dipo table
  first_date = trapping %>% select(Period,date) %>% group_by(Period) %>% summarise_each(funs(min), date)
  dipo_gam = left_join(dipos, first_date, by=c('period'='Period'))
  dipo_gam = rename(dipo_gam, DipoN = x)
  dipo_gam = dipo_gam %>% mutate(month = as.numeric(format(date, "%m")),
                                 Year = as.numeric(format(date, "%Y")),
                                 Time = as.numeric(date) / 1000)
  
  return(dipo_gam)
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