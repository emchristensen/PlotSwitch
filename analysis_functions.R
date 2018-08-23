
#' @title plot GAM prediction
#' @description plot of predictions from GAM model, with data points
#'
#' @param modelPred data frame of predictions from GAM model
#' @param dat data frame of data points used for the GAM
#' @param Palette colours to use
#'
plot_gam_prediction = function(modelPred, dat, Palette) {
  p.plt <- ggplot(modelPred, aes(x = censusdate, y = Fitted)) +
    geom_point(data = dat, mapping = aes(y = n, colour = treatment)) +
    geom_ribbon(aes(ymax = Upper, ymin = Lower, fill = treatment),
                alpha = 0.2) +
    geom_line(aes(colour = treatment)) +
    labs(y = 'Count', x = NULL) +
    theme(legend.position = 'right') +
    scale_colour_manual(name = 'Treatment', values = Palette,
                        breaks=c("CC","EC","XC"),
                        labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
    scale_fill_manual(name = 'Treatment', values = Palette, 
                      breaks=c("CC","EC","XC"),
                      labels=c("long-term\n control", "kangaroo rat\n removal", "rodent\n removal")) +
    geom_vline(xintercept=as.Date('2015-04-10')) 
  return(p.plt)
}


#' @title plot smooth difference
#' @description plot of difference of smoothes from GAM model
#' 
#' @param diffs data frame containing difference of smooth info (output from osmooth_diff())
#' @param Palette colours to use
#'
plot_smooth_diff = function(diffs,Palette) {
  diffPlt <- ggplot(diffs, aes(x = censusdate, y = diff, group = pair, colour = pair)) +
    geom_ribbon(aes(ymax = upper, ymin = lower, colour = NULL, fill = pair), alpha = 0.15) +
    geom_line() +
    #theme(legend.position = 'top') +
    labs(y = 'Control - Treatment', x = NULL) + 
    scale_colour_manual(name = 'Treatment pair', values = Palette,
                        breaks=c('CC-EC','CC-XC'),
                        labels=c('Control - Former \n kangaroo rat removal','Control - Former \n rodent removal')) +
    scale_fill_manual(name = 'Treatment pair', values = Palette,
                      breaks=c('CC-EC','CC-XC'),
                      labels=c('Control - Former \n kangaroo rat removal','Control - Former \n rodent removal')) +
    geom_vline(xintercept=as.Date('2015-04-10')) +
    geom_hline(yintercept=0)
} 


#' @title predict treatment effect
#' @description creates data frame of predictions from GAM model based on the effect of treatment
#'              (excludes plot effect)
#' 
#' @param dat data frame of original data
#' @param np length of prediction you want (how many prediction points)
#' @param MODEL gam model object
#' @param exVars list of variables (from gam) to be excluded in smooth diff
#' 
predict_treat_effect = function(dat, np, MODEL, exVars) {
  # Data to predict at; note the dummy plot - need to provide all variables used to
  # fit the model when predicting
  treatEff <- with(dat,
                   expand.grid(censusdate = seq(min(censusdate), max(censusdate), length = np),
                               treatment  = c('CC','EC','XC'),
                               plot       = 4)) 
  ## create derived variables from the data we want to predict at
  treatEff <- transform(treatEff,
                        oPlot       = ordered(plot),
                        oTreatment  = ordered(treatment, levels = c('CC','EC','XC')),
                        numericdate = as.numeric(censusdate) / 1000)
  
  # actually predict, on link scale so we can get proper CIs, exclude
  treatPred <- as.data.frame(predict(MODEL, treatEff, type = 'link', se.fit = TRUE,
                                     exclude = exVars))
  # bind predictions to data we predicted at
  treatPred <- cbind(treatEff, treatPred)
  # extract inverse of link function from the model
  ilink <- family(MODEL)$linkinv
  # form 95% bayesian credible interval / frequentist across-function confidence interval
  treatPred <- transform(treatPred, Fitted = ilink(fit),
                         Upper = ilink(fit + (2 * se.fit)),
                         Lower = ilink(fit - (2 * se.fit)))
  return(treatPred)
}


#' @title compute differences of smooths
#' @description Compute pairwise differences of smooths when fitted using ordered factors
#' 
#' @param model gam model object
#' @param newdata predicted data (e.g. output from predict_treat_effect)
#' @param smooth_var variable name from 'newdata' to perform smooth upon
#' @param f1 first treatment type for computing difference ('CC','EC', or 'XC')
#' @param f2 second treatment type for computing difference ('CC','EC', or 'XC')
#' @param var variable name from 'newdata' corresponding to ordered treatment
#' @param alpha
#' @param unconditional
#' @param removePara T/F: remove parametric columns not associated with var?
#' @param keepVar
#' 
osmooth_diff <- function(model, newdata, smooth_var, f1, f2, var, alpha = 0.05,
                         unconditional = FALSE, removePara = TRUE, keepVar = TRUE, ...) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix', ...)
  # reference level
  ref_level <- levels(newdata[[var]])[1L]
  ref_smooth <- grepl(paste0("s\\(", smooth_var, "\\)\\.{1}[[:digit:]]+$"), colnames(xp))
  not_smooth <- !grepl('^s\\(', colnames(xp))
  c1 <- ref_smooth | grepl(f1, colnames(xp))
  c2 <- ref_smooth | grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  # difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  # zero out cols of X related to splines for other levels
  X[, !not_smooth][, ! (c1[!not_smooth] | c2[!not_smooth])] <- 0
  if (isTRUE(removePara)) {
    # zero out the parametric cols not associated with `var`,
    # ignore (Intercept), as it is zero anyway
    ind <- grepl('^s\\(', colnames(xp))
    if (isTRUE(keepVar)) {
      ind <- ind | grepl(paste0('^', var), colnames(xp))
    }
    X[, !ind] <- 0
  }
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qnorm(alpha/2, lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  data.frame(pair = paste(f1, f2, sep = '-'),
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr,
             censusdate = newdata[r1,'censusdate'])
}
