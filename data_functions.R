library(dplyr)



#' append dates
#' 
#' @description appends date info based on period
#' 
#' @param df data frame: must have a "period" column
#' 
append_dates = function(df) {
  # Get trapping history
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"
  trapping = read.csv(text=RCurl::getURL(http)) 
  trapping$date = as.Date(paste(trapping$year,trapping$month,trapping$day,sep='-'))
  
  # create data frame with period, date, month, year, time
  period_dates = aggregate(trapping$date,by=list(period=trapping$period),FUN=min)
  names(period_dates) = c('period','date')
  period_dates$month = as.numeric(format(period_dates$date, "%m"))
  period_dates$year = as.numeric(format(period_dates$date, "%Y"))
  period_dates$time = as.numeric(period_dates$date) / 1000
  
  # append to data frame
  df_new = merge(df,period_dates,by='period')
  return(df_new)
}


#' Function to tally captures by species, period, and plot
#' This function was made specifically for working with the 2015 plot switch
#' 
#' @param start_period first period number of data desired; default is 130 (1989)
#' @param incomplete T/F, wheter or not to include incomplete censuses
#' 
#' @return data frame with columns 'plot', 'period', 'species', 'x' (abundance), 'date','month','Year','Time'

rodent_abundance = function(start_period=130,incomplete=F) {
  # filter data: target species, remove early periods, remove incomplete censuses
  rdat_filtered = filter_data(start_period,incomplete)
  
  # aggregate by species, plot, period
  byspecies = aggregate(rdat_filtered$species,by=list(period = rdat_filtered$period, plot = rdat_filtered$plot,species = rdat_filtered$species),FUN=length)
  
  # adding sampling date to the dipo table
  byspecies_date = append_dates(byspecies)
  
  return(byspecies_date)
}


#' @title filter data
#' 
#' @description preliminary step to load and filter data before creating other rodent data metrics
#' 
#' @param
#' 
filter_data = function(start_period,incomplete=F) {
  # load rodent data from repo
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent.csv"
  rdat = read.csv(text=RCurl::getURL(http),as.is=T,na.strings = '')
  
  # Get trapping history in order to find and remove incomplete censuses
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"
  trapping = read.csv(text=RCurl::getURL(http)) 
  plotstrapped = aggregate(trapping$sampled,by=list(period=trapping$period),FUN=sum)
  fullcensus = plotstrapped[plotstrapped$x>=21,] # warning: this may not be ok for other projects, period 457 only trapped 21 plots but the skipped ones aren't relevant to this project
  
  # filter data; remove early trapping periods, nontarget species
  targetsp = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','RO','SF','SH')
  rdat_filtered = dplyr::filter(rdat, period >= start_period, species %in% targetsp)
  
  # if desired, remove incompete trapping periods
  if (incomplete==F) {
    rdat_filtered = filter(rdat_filtered, period %in% fullcensus$period)
  }
  return(rdat_filtered)
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
  byspecies = rodent_abundance(start_period=415,incomplete=F)
  
  # number of dipos per plot
  d = byspecies %>% 
    filter(period>414, species %in% c('DM','DO','DS'))
  dipos = aggregate(d$x, by=list(period=d$period, plot=d$plot, date=d$date, month=d$month, year=d$year, time=d$time), FUN=sum)
  # data frame of all plots in all periods
  allplotsperiod = expand.grid(period=unique(dipos$period), plot=unique(dipos$plot))
  #attach date columns according to period
  allplotsperiod = append_dates(allplotsperiod)
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
  dipos_gam = rename(dipos,DipoN=x)
  # put data in chronological order
  dipos_gam = dipos_gam[order(dipos_gam$period),]
  
  return(dipos_gam)
}

#' @title Make Data
#'
#' @description make a data frame for running GAM models
#' 
#' @param species species desired; either "All," "Granivore," "Dipos," "SmGran," "SmH" (small heteromyid), or "SmM" (small murid)
#' @param start_period first period number of data desired; default is 415 (2013)
#' 
#' @return data frame with columns:
#'            plot
#'            period
#'            treatment (two letters representing before/after switch: C = control, E = krat exclosure, X = total rodent removal)
#'            n (number of small heteromyids)
#'            date 
#'            month
#'            Year
#'            Time

make_data= function(species='All',start_period=415) {
  # Create species abundances by plot by period
  dat = rodent_abundance(start_period=start_period,incomplete=F)
  # select species of interest
  if (species=='All') {targetsp = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','RO','SF','SH')}
  if (species=='Granivore') {targetsp = c('BA','DM','DO','DS','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')}
  if (species=='SmGran') {targetsp = c('BA','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')}
  if (species=='SmH') {targetsp = c('PB','PF','PH','PI','PP')}
  if (species=='SmM') {targetsp = c('BA','PE','PL','PM','RF','RM','RO')}
  if (species=='Dipos') {targetsp = c('DM','DO','DS')}
  
  target_dat = filter(dat,species %in% targetsp)
  total = aggregate(target_dat$x,by=list(period=target_dat$period,plot=target_dat$plot,date=target_dat$date,month=target_dat$month,
                                         year=target_dat$year,time=target_dat$time),FUN=sum)

  # data frame of all plots in all periods
  allplotsperiod = expand.grid(period=unique(dat$period), plot=unique(dat$plot))
  #attach date columns according to period
  allplotsperiod = append_dates(allplotsperiod)
  # attach treatment column according to plot number
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  allplotsperiod = merge(allplotsperiod,treatment)
  # merge capture data with data frame of all plots and all periods, and fill in empty data with zeros
  total = merge(allplotsperiod,total,all=T)
  total$x[is.na(total$x)] = 0
  # change column name
  total_gam = rename(total,n=x)
  # put data in chronological order
  total_gam = total_gam[order(total_gam$period),]
  
  return(total_gam)
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
  CE = data %>% filter(treatment == "CE") %>% arrange(date)
  EE = data %>% filter(treatment == "EE") %>% arrange(date)
  return(list(CC,EC,XC,CE,EE))
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
                       by=list(period=rdat$period,plot=rdat$plot),
                       FUN=unique)
  for (n in 1:length(richness$period)) {
    richness$nsp[n] = length(unlist(richness$x[n]))
  }
  
  # make sure there are zeros where plot was trapped but nothing caught; add column for treatment type
  allplotsperiod = expand.grid(period=unique(richness$period), plot=unique(richness$plot))
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  allplotsperiod = merge(allplotsperiod,treatment)
  sprich = merge(allplotsperiod,richness,all=T)
  sprich$nsp[is.na(sprich$nsp)] <- 0
  sprich = append_dates(sprich)
  # put data in chronological order
  sprich = sprich[order(sprich$period),]
  
  return(dplyr::select(sprich,period,plot,nsp,treatment,date,month,year,time))
}

#' @title get mass data
#' 
#' @description get total or average rodent mass by period and plot
#'
#' @param start_period first period number of data desired; default is 130 (1989)
#' @param avg T/F: if true, computes average wgt per plot, if F computes total wgt
#' @param metE T/F: if T converts mass to metabolic energy, if F returns mass data
#' 
#' @example avg_mass = get_mass(start_period = 415,avg=T)
#'
#'
get_mass = function(start_period=130, avg, metE) {
  # load rodent data
  rdat_filtered = filter_data(start_period,incomplete=F) %>% select(period,plot,species,wgt)
  
  # fill in unweighed rodents with avg mass for that species
  avgmass = aggregate(rdat_filtered$wgt, by = list(species = rdat_filtered$species),FUN=mean,na.rm=T)
  rdat_filtered = merge(rdat_filtered,avgmass)
  for (n in seq(length(rdat_filtered$species))) {
    if (is.na(rdat_filtered$wgt[n])) {
      rdat_filtered$wgt[n] = rdat_filtered$x[n]
    }
  }
  
  # convert mass to metabolic energy if desired
  if (metE == T) {
    rdat_filtered$wgt = 5.69 * rdat_filtered$wgt^0.75
  }
  
  # get either total or average mass per plot and period
  if (avg == T) {
    mass = aggregate(rdat_filtered$wgt,by=list(period=rdat_filtered$period,plot=rdat_filtered$plot),FUN=mean)
  } else {
    mass = aggregate(rdat_filtered$wgt,by=list(period=rdat_filtered$period,plot=rdat_filtered$plot),FUN=sum)
  }
  
  # change column name
  mass = rename(mass,n=x)
  
  # data frame of all plots in all periods
  allplotsperiod = expand.grid(period=unique(rdat_filtered$period), plot=unique(rdat_filtered$plot))
  #attach date columns according to period
  allplotsperiod = append_dates(allplotsperiod)
  # attach treatment column according to plot number
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  allplotsperiod = merge(allplotsperiod,treatment)
  # merge capture data with data frame of all plots and all periods, and fill in empty data with zeros
  mass_data = merge(allplotsperiod,mass,all=T)
  mass_data$n[is.na(mass_data$n)] = 0
  # put data in chronological order
  mass_data = mass_data[order(mass_data$period),]
  return(mass_data)
}

#' @title get species evenness
#' 
#' @description calculate species evenness per plot/period
#' 
#' @param start_period first period number of data desired
#' 
get_evenness = function(start_period=415) {
  # load rodent data
  rdat_filtered = filter_data(start_period,incomplete=F) %>% select('period','plot','species')
  
  # put data in tabular form
  rdat_filtered$periodplot = paste(rdat_filtered$period,rdat_filtered$plot)
  plt_table = table(rdat_filtered$periodplot,rdat_filtered$species)

  # calculate Shannon diversity and evenness
  H = vegan::diversity(plt_table)
  even = H/log(vegan::specnumber(plt_table))
  
  evenness = data.frame(periodplot = names(even), n = unname(even))
  evenness$periodplot = as.character(evenness$periodplot)
  evenness$period = rep(NA)
  evenness$plot = rep(NA)
  for (ind in seq(length(evenness$periodplot))) {
    evenness$period[ind] = strsplit(evenness$periodplot,' ')[[ind]][1]
    evenness$plot[ind] = strsplit(evenness$periodplot,' ')[[ind]][2]
  }
  # append date columns
  evenness_dat = evenness %>% select(-periodplot) %>% append_dates()
  # attach treatment column according to plot number
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  evenness_dat = merge(evenness_dat,treatment)
  # put rows in chronological order
  evenness_dat = evenness_dat[order(evenness_dat$period),]
  
  return(evenness_dat)
}
