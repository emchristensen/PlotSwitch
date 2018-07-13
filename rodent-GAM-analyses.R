library(dplyr)
library(mgcv)
library(ggplot2)

source('figure_functions.R')

# ==========================================================================================
# Number of dipodomys
rodent <- read.csv('Dipo_counts.csv')
rodent$censusdate <-as.Date(rodent$censusdate)

# create variables needed for GAM
rodent <- dplyr::mutate(rodent,
                 oTreatment = ordered(treatment, levels = c('CC','EC','XC')),
                 oPlot      = ordered(plot),
                 plot       = factor(plot))

# control
ctrl <- gam.control(nthreads = 4)

# GAM model --- includes plot-specific smooth differences
dipo.gam <- gam(n ~ oPlot + oTreatment + s(numericdate, k = 20) +
                  s(numericdate, by = oTreatment, k = 15) +
                  s(numericdate, by = oPlot),
                data = rodent, method = 'REML', family = poisson, select = TRUE, control = ctrl)


# Look at the treatment effect smooths on count scale. 
# Requires us to exclude plot effects
treatPred <- predict_treat_effect(rodent, np = 500, MODEL=dipo.gam)

# extract inverse of link function from the model
ilink <- family(dipo.gam)$linkinv

# plot GAM fit and data
p.plt = plot_gam_prediction(treatPred,rodent)
p.plt

#ggsave('estimated-treatment-effects.pdf', p.plt)


# Compute pairwise treatment diffs if we leave *in* the parametric Treatment terms
d1 <- osmooth_diff(MODEL, treatEff, "numericdate", "CC", "EC", var = "oTreatment", removePara = FALSE)
d2 <- osmooth_diff(MODEL, treatEff, "numericdate", "CC", "XC", var = "oTreatment", removePara = FALSE)

diffs <- rbind(d1, d2)

## plot difference
diffPlt = plot_smooth_diff(diffs)
diffPlt
