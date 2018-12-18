library(dplyr)
library(mgcv)
library(ggplot2)
library(cowplot)

source('FinalAnalysis/analysis_functions.R')
theme_set(theme_bw())
cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")

# ==========================================================================================
# Number of dipodomys
dipo <- read.csv('Data/Dipo_counts.csv')
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
# terms to exclude; must be named exactly as printed in `summary(model)` output
exVars.d <- c('oPlot', paste0('s(numericdate):oPlot', c(5,6,7,11,13,14,17,18,24)))
treatPred.dipo <- predict_treat_effect(dipo, np = 500, MODEL=dipo.gam, exVars.d)

# plot GAM fit and data
d.plt <- plot_gam_prediction(treatPred.dipo, dipo, Palette=cbPalette[1:3], ylab='Count')
d.plt

#ggsave('Figures/dipo-treatment-effects.png', d.plt,width=6,height=2.5)

# Compute pairwise treatment diffs if we leave *in* the parametric Treatment terms
d1 <- osmooth_diff(dipo.gam, treatPred.dipo, "numericdate", "CC", "EC", var = "oTreatment", removePara = FALSE)
d2 <- osmooth_diff(dipo.gam, treatPred.dipo, "numericdate", "CC", "XC", var = "oTreatment", removePara = FALSE)
diffs.dipo <- rbind(d1, d2)

## difference of smooths
diffPlt <- plot_smooth_diff(diffs.dipo, Palette=cbPalette[2:3])
diffPlt
#ggsave('dipo-difference.png', diffPlt,width=6,height=2.5)

## Cowplot grid
dipo_plot = plot_grid(d.plt, diffPlt, labels = "AUTO", ncol = 1, align = 'v')
dipo_plot
ggsave('Figures/dipo-gam-plots.png', dipo_plot, width=4, height = 4.2, dpi=300)

# =========================================================================================
# number of small granivores
smgran <- read.csv('Data/SmallGranivores.csv')
smgran$censusdate <- as.Date(smgran$censusdate)

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
# terms to exclude; must be named exactly as printed in `summary(model)` output
exVars.sg <- c('oPlot', paste0('s(numericdate):oPlot', c(5,6,7,11,13,14,17,18,24)))
treatPred.sg <- predict_treat_effect(smgran, np = 500, MODEL=smgran.gam, exVars.sg)

sg.plt <- plot_gam_prediction(treatPred.sg, smgran, Palette=cbPalette[1:3], ylab='Count')
sg.plt
#ggsave('smallgran-treatment-effects.png', sg.plt,width=6,height=2.5)

# difference of smooths
d1 <- osmooth_diff(smgran.gam, treatPred.sg, "numericdate", "CC", "EC", var = "oTreatment", removePara = FALSE)
d2 <- osmooth_diff(smgran.gam, treatPred.sg, "numericdate", "CC", "XC", var = "oTreatment", removePara = FALSE)
diffs.sg <- rbind(d1,d2)
sg.diffPlt <- plot_smooth_diff(diffs.sg,Palette=cbPalette[2:3])
sg.diffPlt
#ggsave('smallgran-difference.png', sg.diffPlt,width=6,height=2.5)

## Cowplot grid
sg_plot = plot_grid(sg.plt, sg.diffPlt, labels = "AUTO", ncol = 1, align = 'v')
sg_plot
ggsave('Figures/smallgran-gam-plots.png', sg_plot, width=4, height = 4.2, dpi=300)

# ========================================================================================
# Total rodent energy
energy <- read.csv('Data/TotalCommunityEnergy.csv')
energy$censusdate <- as.Date(energy$censusdate)

energy <- mutate(energy,
                 oTreatment = ordered(treatment, levels = c('CC','EC','XC')),
                 oPlot      = ordered(plot),
                 plot       = factor(plot))

# GAM model -- includes plot smooths
energy.gam <- gam(n ~ oPlot + oTreatment + s(numericdate, k = 20) +
                    s(numericdate, by = oTreatment, k = 15) +
                    s(numericdate, by = oPlot),
                  data = energy, method = 'REML', family = tw, select = TRUE)

# plot treatment effects (exclude plot smooths)
exVars.energy <- c('oPlot', paste0('s(numericdate):oPlot', c(5,6,7,11,13,14,17,18,24)))
treatPred.energy <- predict_treat_effect(energy, np = 500, MODEL=energy.gam, exVars.energy)

energy.plt <- plot_gam_prediction(treatPred.energy, energy, Palette = cbPalette[1:3], ylab='Metabolic Flux \n(kJ)')
energy.plt
#ggsave('energy-treatment-effects.png', energy.plt,width=6,height=2.5)

# difference of smooths
d1 <- osmooth_diff(energy.gam, treatPred.energy, "numericdate", "CC", "EC", var = "oTreatment", removePara = FALSE)
d2 <- osmooth_diff(energy.gam, treatPred.energy, "numericdate", "CC", "XC", var = "oTreatment", removePara = FALSE)
diffs.energy <- rbind(d1,d2)
energy.diffPlt <- plot_smooth_diff(diffs.energy,Palette=cbPalette[2:3])
energy.diffPlt
#ggsave('energy-difference.png', energy.diffPlt,width=6,height=2.5)

## Cowplot grid
energy_plot = plot_grid(energy.plt, energy.diffPlt, labels = "AUTO", ncol = 1, align = 'v')
energy_plot
ggsave('Figures/energy-gam-plots.png', energy_plot, width=4, height = 4.2, dpi=300)

# # ==========================================================================================
# # Species richness
# sprich <- read.csv('Data/SpeciesRichness.csv')
# sprich$censusdate <- as.Date(sprich$censusdate)
# 
# sprich <- mutate(sprich,
#                  oTreatment = ordered(treatment, levels = c('CC','EC','XC')),
#                  oPlot      = ordered(plot),
#                  plot       = factor(plot))
# 
# # model 1 --- this has an intercept for plot but plots follow respective treatment smooth
# sprich.gam <- gam(n ~ oTreatment + s(numericdate, k = 20) +
#                     s(numericdate, by = oTreatment, k = 15) +
#                     s(plot, bs = "re"),
#                   data = sprich, method = 'REML', family = poisson, select = TRUE, control = gam.control(nthreads = 4))
# 
# # plot treatment effects
# # terms to exclude
# exVars.rich <- 's(plot)'
# treatPred.rich <- predict_treat_effect(sprich, np = 500, MODEL=sprich.gam, exVars.rich)
# 
# rich.plt <- plot_gam_prediction(treatPred.rich, sprich, Palette = cbPalette[1:3], ylab='# species')
# rich.plt
# #ggsave('richness-treatment-effects.png', rich.plt,width=6,height=2.5)
# 
# # difference of smooths
# d1 <- osmooth_diff(sprich.gam, treatPred.rich, "numericdate", "CC", "EC", var = "oTreatment", removePara = FALSE)
# d2 <- osmooth_diff(sprich.gam, treatPred.rich, "numericdate", "CC", "XC", var = "oTreatment", removePara = FALSE)
# diffs.rich <- rbind(d1,d2)
# rich.diffPlt <- plot_smooth_diff(diffs.rich,Palette=cbPalette[2:3])
# rich.diffPlt
# #ggsave('richness-difference.png', rich.diffPlt,width=6,height=2.5)
# 
# ## Cowplot grid
# rich_plot = plot_grid(rich.plt, rich.diffPlt, labels = "AUTO", ncol = 1, align = 'v')
# rich_plot
# #ggsave('Figures/sprich-gam-plots.png', rich_plot, width=7, height = 5)
# 
