library(dplyr)

source('rodent_abundance_by_period_and_plot.r')

make_prediction = function(data, model, numpredictions=50){
  
  # Make data to plot the trend line on the data
  avg = mean(data$DipoN)
  want <- seq(1, nrow(data), length.out = 50)
  pdat <- with(data, data.frame(Time = Time[want], date = date[want], 
                                month = month[want]))
  p  <- predict(model$gam,  newdata = pdat, type = "terms", se.fit = TRUE)
  pdat <- transform(pdat, p  = p$fit[,2], p_raw=p$fit[,2] +avg, se  = p$se.fit[,2])
  return(pdat)
}

make_dipo_data = function(){
  # Create species abundances by plot by period
  
  byspecies = rodent_abundance('Granivore')
  treatment = data.frame(treatment = c('CX','CE','EE','CC',
                                       'XC','EC','XC','CE',
                                       'CX','XX','CC','CX',
                                       'EC','CC','EE','XX',
                                       'CC','EC','EE','EE',
                                       'EE','CE','XX','XC'),plot=seq(1,24))
  http = "https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv"
  trapping = read.csv(text=getURL(http)) 
  trapping$date = as.Date(paste(trapping$Year,trapping$Month,trapping$Day,sep='-'))
  plotstrapped = aggregate(trapping$Sampled,by=list(period=trapping$Period),FUN=sum)
  fullcensus = plotstrapped[plotstrapped$x>20,]
  
  byspecies$type = paste(byspecies$type_before,byspecies$type_after,sep='')
  
  rdat = byspecies %>% filter(period>414,type %in% c('EC','CC','XC'),period %in% fullcensus$period) %>% select(plot,period,species,x,type)
  
  # number of dipos per plot
  
  d = filter(rdat,species %in% c('DM','DO','DS'))
  dipos = aggregate(d$x,by=list(period=d$period,treatment=d$type,plot=d$plot),FUN=sum)
  allplotsperiod = expand.grid(period=unique(dipos$period),plot=unique(dipos$plot))
  allplotsperiod = merge(allplotsperiod,treatment)
  
  dipos = merge(allplotsperiod,dipos,all=T)
  dipos$x[is.na(dipos$x)] = 0
  
  # adding sampling date to the dipo table
  
  first_date = trapping %>% select(Period,date) %>% group_by(Period) %>% summarise_each(funs(min), date)
  first_date = rename(first_date,period=Period)
  
  dipo_gam = left_join(dipos, first_date)
  dipo_gam = rename(dipo_gam, DipoN = x)
  dipo_gam = dipo_gam %>% mutate(month = as.numeric(format(date, "%m")),
                     Year = as.numeric(format(date, "%Y")),
                     Time = as.numeric(date) / 1000)
  return(dipo_gam)
}

gam_diagnostics = function(model, title){
  
  layout(matrix(1:2, ncol = 2))
  plot(model$gam, scale = 0, main = title)
  layout(1)
  
  layout(matrix(1:2, ncol = 2))
  acf(resid(model$lme), lag.max = 36, main = "ACF")
  pacf(resid(model$lme), lag.max = 36, main = "pACF")
  layout(1)
  return(summary(model$gam))
}