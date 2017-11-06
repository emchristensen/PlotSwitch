library(dplyr)
library(portalr)

#################################################
# TO DO:
# 1) In the original code for creating GAM data, period 457 was added back into the
#    data. It is an incomplete census, but all the plots needed for this analysis
#    were sampled. In refactoring the code to work with portalr, SKME did not
#    add back in 457 only because doing so would be a pain. Adding 457 still 
#    needs to be done. Below is the orginal line of code on this issue. It won't
#    work not because of the way portalr is being used:
#    fullcensus = plotstrapped[plotstrapped$x>=21,] 
#    warning: this may not be ok for other projects, period 457 only trapped 21 
#    plots but the skipped ones aren't relevant to this project

# test code: 
# dat = get_data(startdate = "2013-03-11")
# dat = make_data(species="SmM")
#################################################

#' @title Get Rodent Data
#' 
#' @description make a data frame for running GAM models
#' 
#' @param startdate Census date of the period code used to filter the data
#' 
#' @return data frame with columns:
#'            plot
#'            censusdate
#'            species
#'            abundance
#'            numericdate
#'            treatment (two letters representing before/after switch: C = control, E = krat exclosure, X = total rodent removal)


get_data = function(startdate){
  data = portalr::abundance(path='repo', level = 'Plot', type='Granivores',
                            length="All", unknowns=TRUE, incomplete=FALSE,
                            shape="flat", time='date')
  data$numericdate = as.numeric(data$censusdate) / 1000
  rdat_filtered = dplyr::filter(data, censusdate >= startdate)
  rdat_filtered = add_treatment(rdat_filtered)
}

#' @title Make Data
#'
#' @description make a data frame of rodent abundances for GAM models
#' 
#' @param species species desired; either "All," "Granivore," "Dipos," "SmGran," "SmH" (small heteromyid), or "SmM" (small murid)
#' @param start_period first period number of data desired; default is 415 (2013)
#' 
#' @return data frame with columns:
#'            plot
#'            censusdate (date of data collection formatted YYYY-MM-DD)
#'            treatment (two letters representing before/after switch: C = control, E = krat exclosure, X = total rodent removal)
#'            n (number of individuals)
#'            Time (numeric version of censusdate for GAM analyses)

make_data= function(species='All',start_date = "2013-03-11") {
  # Create species abundances by plot by period
  dat = get_data(startdate=start_date)
  # select species of interest
  if (species=='All') {targetsp = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','RO','SF','SH')}
  if (species=='Granivore') {targetsp = c('BA','DM','DO','DS','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')}
  if (species=='SmGran') {targetsp = c('BA','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO')}
  if (species=='SmH') {targetsp = c('PB','PF','PH','PI','PP')}
  if (species=='SmM') {targetsp = c('BA','PE','PL','PM','RF','RM','RO')}
  if (species=='Dipos') {targetsp = c('DM','DO','DS')}
  
  target_dat = dplyr::filter(dat,species %in% targetsp)
  total = aggregate(target_dat$abundance,
                    by=list(plot=target_dat$plot,censusdate=target_dat$censusdate,
                            Time=target_dat$numericdate, 
                            treatment=target_dat$treatment),
                    FUN=sum)
  
  total$x[is.na(total$x)] = 0
  # change column name
  total_gam = rename(total,n=x)
  # put data in chronological order
  total_gam = total_gam[order(total_gam$Time),]
  trt_files = trt_data(total_gam)
  
  return(trt_files)
}

#' @title Add treatment codes
#' 
#' @description adds treatment codes for each plot in the rodent dataframe
#' 
#' @param data dataframe must have a column named 'plot'
#' 
#' @return dataframe with columns:
#'            plot
#'            censusdate
#'            species
#'            abundance
#'            numericdate
#'            treatment (two letters representing before/after switch: C = control, E = krat exclosure, X = total rodent removal)
add_treatment = function(data){
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  data = merge(data,treatment)
  data = data %>% dplyr::filter(treatment %in% c('CC','XC','EC'))
  return(data)
}


#' @title Separate Treatment Data
#' 
#' @description takes data with muyltiple treatment and separates into distinct files for analysis
#' 
#' @param data dataframe output from make_data
#' 
#' @return a list containing three data frames, containing dipo data from 3 treatment types

trt_data = function(data){
  data$plot = as.factor(data$plot)
  CC = dplyr::filter(data, treatment == "CC") %>% dplyr::arrange(Time)
  EC = dplyr::filter(data, treatment == "EC") %>% dplyr::arrange(Time)
  XC = dplyr::filter(data, treatment == "XC") %>% dplyr::arrange(Time)
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
  
  return(dplyr::select(sprich,period,plot,nsp,treatment,date,month,year,Time))
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

make_dipo_data = function(){
  # Create species abundances by plot by period
  byspecies = rodent_abundance(start_period=415,incomplete=F)
  
  # number of dipos per plot
  d = byspecies %>% 
    filter(period>414, species %in% c('DM','DO','DS'))
  dipos = aggregate(d$x, by=list(period=d$period, plot=d$plot, date=d$date, month=d$month, year=d$year, Time=d$Time), FUN=sum)
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