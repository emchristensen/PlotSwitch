
#' @description 
#' 
#' @param data
#' @param model
#' @param numpredictions
#' 
#' 
make_prediction_gamm = function(data, model, numpredictions=50){
  
  # Make data to plot the trend line on the data
  avg = mean(data$DipoN)
  times = select(data,date,month,Time) %>% unique()
  want <- seq(1, nrow(times), length.out = 50)
  pdat <- with(times, data.frame(Time = Time[want], date = date[want], 
                                 month = month[want]))
  p  <- predict(model$gam,  newdata = pdat, type = "terms", se.fit = TRUE)
  pdat <- transform(pdat, p  = p$fit[,2], p_raw=p$fit[,2] +avg, se  = p$se.fit[,2])
  pdat = pdat %>% mutate(lower = p_raw - 1.96*se, upper = p_raw + 1.96*se)
  return(pdat)
}

gamm_diagnostics = function(model, title){
  
  layout(matrix(1:2, ncol = 2))
  plot(model$gam, scale = 0, main = title)
  layout(1)
  
  layout(matrix(1:2, ncol = 2))
  acf(resid(model$lme), lag.max = 36, main = "ACF")
  pacf(resid(model$lme), lag.max = 36, main = "pACF")
  layout(1)
  return(summary(model$gam))
}