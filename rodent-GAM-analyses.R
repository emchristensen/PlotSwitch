library(dplyr)
library(mgcv)
library(ggplot2)

source('analysis_functions.R')
theme_set(theme_bw())

# ==========================================================================================
# Number of dipodomys
dipo <- read.csv('Dipo_counts.csv')
dipo$censusdate <-as.Date(dipo$censusdate)

# create variables needed for GAM
dipo <- dplyr::mutate(dipo,
                 oTreatment = ordered(treatment, levels = c('CC','EC','XC')),
                 oPlot      = ordered(plot),
                 plot       = factor(plot))

# GAM model --- includes plot-specific smooths
dipo.gam <- gam(n ~ oPlot + oTreatment + s(numericdate, k = 20) +
                  s(numericdate, by = oTreatment, k = 15) +
                  s(numericdate, by = oPlot),
                data = dipo, method = 'REML', family = poisson, select = TRUE, control = gam.control(nthreads = 4))

# Look at the treatment effect smooths on count scale. 
# Requires us to exclude plot effects
treatPred.dipo <- predict_treat_effect(dipo, np = 500, MODEL=dipo.gam)

# plot GAM fit and data
d.plt = plot_gam_prediction(treatPred.dipo,dipo)
d.plt

#ggsave('estimated-treatment-effects.pdf', p.plt)

# Compute pairwise treatment diffs if we leave *in* the parametric Treatment terms
d1 <- osmooth_diff(dipo.gam, treatPred.dipo, "numericdate", "CC", "EC", var = "oTreatment", removePara = FALSE)
d2 <- osmooth_diff(dipo.gam, treatPred.dipo, "numericdate", "CC", "XC", var = "oTreatment", removePara = FALSE)
diffs.dipo <- rbind(d1, d2)

## difference of smooths
diffPlt = plot_smooth_diff(diffs.dipo)
diffPlt

# =========================================================================================
# number of small granivores
smgran <- read.csv('SmallGranivores.csv')
smgran$censusdate = as.Date(smgran$censusdate)

smgran <- mutate(smgran,
                 oTreatment = ordered(treatment, levels = c('CC','EC','XC')),
                 oPlot      = ordered(plot),
                 plot       = factor(plot))

# GAM model - plot and treatment smooths
smgran.gam <- gam(n ~ oPlot + oTreatment + s(numericdate, k = 20) +
                    s(numericdate, by = oTreatment, k = 15) +
                    s(numericdate, by = oPlot),
                  data = smgran, method = 'REML', family = poisson, select = TRUE, control = gam.control(nthreads = 4))

# plot treatment effects
treatPred.sg = predict_treat_effect(smgran, np = 500, MODEL=smgran.gam)

sg.plt = plot_gam_prediction(treatPred.sg, smgran)
sg.plt

# difference of smooths
d1 = osmooth_diff(smgran.gam, treatPred.sg, "numericdate", "CC", "EC", var = "oTreatment", removePara = FALSE)
d2 = osmooth_diff(smgran.gam, treatPred.sg, "numericdate", "CC", "XC", var = "oTreatment", removePara = FALSE)
diffs.sg = rbind(d1,d2)
sg.diffPlt = plot_smooth_diff(diffs.sg)
sg.diffPlt
