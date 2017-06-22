library(dplyr)



make_prediction_gam = function(data, model, numpredictions=50){
  
  # Make data to plot the trend line on the data
  avg = mean(data$DipoN)
  times = select(data,date,month,Time) %>% unique()
  want <- seq(1, nrow(times), length.out = 50)
  pdat <- with(times, data.frame(Time = Time[want], date = date[want], 
                                 month = month[want]))
  p  <- predict(model,  newdata = pdat, type = "terms", se.fit = TRUE)
  pdat <- transform(pdat, p  = p$fit[,2], p_raw=p$fit[,2] +avg, se  = p$se.fit[,2])
  pdat = pdat %>% mutate(lower = p_raw - 1.96*se, upper = p_raw + 1.96*se)
  return(pdat)
}


gam_diagnostics = function(model, title){
  
  layout(matrix(1:2, ncol = 2))
  plot(model, scale = 0, main = title)
  layout(1)
  
  layout(matrix(1:2, ncol = 2))
  acf(resid(model), lag.max = 36, main = "ACF")
  pacf(resid(model), lag.max = 36, main = "pACF")
  layout(1)
  return(summary(model))
}
  
plot_singleGAM = function(data, title, ylab, treatment){
  if (treatment == 'CC'){
    lincol= 'blue'
    title = paste("Control -", title, sep= " ")
  } else if (treatment == 'XC'){
    lincol='green'
    title = paste("Krat Exclosure -", title, sep= " ")
  } else {
    lincol='red'
    title = paste("Removals -", title, sep= " ")
    }
  transition = as.Date("2015-03-15", format="%Y-%m-%d")
  ggplot(aes(x=date, y=p_raw), data = data) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill='gray90') +
    geom_line(color=lincol) +
    geom_vline(xintercept =  as.numeric(transition)) +
    ggtitle(title) +
    xlab("Date") + ylab("Dipodomys abundance per plot") +
    theme_classic()
}
