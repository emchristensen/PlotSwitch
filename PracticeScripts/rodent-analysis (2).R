## Analyse the full rodent data

## packages
library('mgcv')
library('ggplot2')
library('readr')
## install.packages('countreg', repos = 'http://R-Forge.R-project.org')
library('countreg')
library('dplyr')
library('tidyr')
library('cowplot')
theme_set(theme_bw())

## load data
rodent <- read_csv('Dipo_counts_2017.csv')

## create variables needed
rodent <- mutate(rodent,
                 oTreatment = ordered(treatment, levels = c('CC','EC','XC')),
                 oPlot      = ordered(plot),
                 plot       = factor(plot))

## model 1 --- this has an intercept for plot but plots follow respective treatment smooth
## truend select = TRUE on here to provide some extra regularisation as we need it for the
## more complex model...
m1 <- gam(dipos ~ oTreatment + s(numericdate, k = 20) +
              s(numericdate, by = oTreatment, k = 15) +
              s(plot, bs = "re"),
          data = rodent, method = 'REML', family = poisson, select = TRUE)

gam.check(m1)          # not too bad, resids are iffy, but low means
summary(m1)            # looks like we don't need plot-level intercepts
plot(m1, pages = 1, shade = TRUE, scale = 0)

## model 2 --- this is similar to above but now we add plot-specific smooth differences
## need `select` here as there is a little issue with identifiability as plot is nested
## in treatment. SE explodes for one treatment at very end without regularisation
m2 <- gam(dipos ~ oPlot + oTreatment + s(numericdate, k = 20) +
              s(numericdate, by = oTreatment, k = 15) +
              s(numericdate, by = oPlot),
          data = rodent, method = 'REML', family = poisson, select = TRUE)

gam.check(m2, rep = 200) # not too bad, resids are iffy, but low mean counts
summary(m2)       # most plot-specific differences are shrunk to 0 EDF, except 14, 17, 13, 5
plot(m2, pages = 1, shade = TRUE, scale = 0)

AIC(m1, m2)                             # favours plot-specific smooths

## rootograms help diagnose model fits
## compare observed count distribution with the expected counts from the
## conditional distribution of the response --- here Poisson
root.m1 <- rootogram(m1, style = 'hanging', plot = FALSE)
root.m2 <- rootogram(m2, style = 'hanging', plot = FALSE)

## these aren't too bad; over-predicting at low counts (1s, 2s) and over-predicting
## by similar amounts at intermediate counts (4s, 5s).
autoplot(root.m1)
autoplot(root.m2)

## I thought maybe there might be some modest zero-inflation but a quick check with
## a family = ziP fit made no difference. If any zero-inflation is present, it is
## likely to be time-dependent so that needs a more complex fit that I dind't explore
## here.

## See if we can do better on the counts with a NegBin response
## model 3: as model 1 but with negbin
m3 <- gam(dipos ~ oTreatment + s(numericdate, k = 20) +
              s(numericdate, by = oTreatment, k = 15) +
              s(plot, bs = "re"),
          data = rodent, method = 'REML', family = nb, select = TRUE)

gam.check(m3)
summary(m3) # note massive theta == tiny over dispersion rel. to Possion
plot(m3, pages = 1, shade = TRUE, scale = 0)

## model 4: as model 2 but with negbin
m4 <- gam(dipos ~ oPlot + oTreatment + s(numericdate, k = 20) +
              s(numericdate, by = oTreatment, k = 15) +
              s(numericdate, by = oPlot),
          data = rodent, method = 'REML', family = nb, select = TRUE)

gam.check(m4)
summary(m4) # note massive theta == tiny over dispersion rel. to Possion
plot(m4, pages = 1, shade = TRUE, scale = 0)

AIC(m3, m4) # as these are just Poisson fits, m4 beats m2 as we saw earlier

## rootograms
root.m3 <- rootogram(m3, style = 'hanging', plot = FALSE)
root.m4 <- rootogram(m4, style = 'hanging', plot = FALSE)

## no change
autoplot(root.m3)
autoplot(root.m4)

AIC(m1, m2, m3, m4)   # Poisson's give better fits; theta is redundant

## Look at full model fits and compare with observations
## Data to predict at
np <- 500
newd <- with(rodent,
             expand.grid(plot       = unique(plot),
                         censusdate = seq(min(censusdate), max(censusdate), length = np)))
## need a lookup table to get treatment for each plot
lookup <- unique(rodent[c('plot', 'treatment')])
## merge prediction data with treatment lookup
newd <- left_join(newd, lookup)
## add some derived variables needed for model
newd <- transform(newd,
                  oPlot       = ordered(plot),
                  oTreatment  = ordered(treatment, levels = c('CC','EC','XC')),
                  numericdate = as.numeric(censusdate) / 1000)

## predict from all the models...
p1 <- predict(m1, newd, type = 'response')
p2 <- predict(m2, newd, type = 'response')
p3 <- predict(m3, newd, type = 'response')
p4 <- predict(m4, newd, type = 'response')
## ...and bind on to the prediction data
preds <- cbind(newd, p1, p2, p3, p4)

## need to stack this tidily for ggplot
pdata <- gather(preds, Model, Fitted, -(1:6))

## Add variables for model family and model form;
## Type 1 is model with only a plot-level intercept
## Type 2 is model with plot level intercept plus plot-level difference smooths
pdata <- mutate(pdata,
                Family = case_when(Model %in% c('p1', 'p2') ~ 'Poisson',
                                   Model %in% c('p3', 'p4') ~ 'Negative Binomial'),
                ModelForm = case_when(Model %in% c('p1', 'p3') ~ 'Type 1',
                                      Model %in% c('p2', 'p4') ~ 'Type 2'))

## visualise
ggplot(pdata, aes(x = censusdate, y = Fitted, colour = ModelForm)) +
    geom_point(aes(x = censusdate, y = dipos), data = rodent, inherit.aes = FALSE) +
    geom_line(size = 1) +
    facet_grid(plot ~ Family) +
    labs(y = 'Count', x = NULL) +
    theme(legend.position = 'top') +
    scale_colour_brewer(type = 'qual', palette = 'Dark2')

## So plot-level smooth differences (i.e. plot-specific trends superimposed on
## treatment-specific smooths make modest differences, at some sites only.
## AIC suggests the more flexible model (Type 2) has better fit, and summary()
## tests indicate that at least for some plots the effect is different from a 0
## effect.

## Also, there's no difference between NB and Poisson fits. You can see this in the
## summary() output for any of the `family = nb` models. Look for the value of theta
## at the top of the summary output --- it's huge, which translates, given the
## parametrisation used into a tiny NegBin variance, which means overdispersion
## "effect" is tiny ---> poisson is fine.

## Look at the treatment effect smooths on count scale. Requires us to exclude the
## plot-level effects.

## Data to predict at; note the dummy plot - need to provide all variables used to
## fit the model when predicting
treatEff <- with(rodent,
                 expand.grid(censusdate = seq(min(censusdate), max(censusdate), length = np),
                             treatment  = c('CC','EC','XC'),
                             plot       = 4))
## create derived variables from the data we want to predict at
treatEff <- transform(treatEff,
                      oPlot       = ordered(plot),
                      oTreatment  = ordered(treatment, levels = c('CC','EC','XC')),
                      numericdate = as.numeric(censusdate) / 1000)

## terms to exclude; must be named exactly as printed in `summary(model)` output
exVars.m1 <- 's(plot)'                  # this is for m1
exVars.m2 <- c('oPlot', paste0('s(numericdate):oPlot', c(5,6,7,11,13,14,17,18,24))) # m2

## Choose which model; Poisson response is fine, and
## model with smooth plot-level differences fits slightly better == m2
MODEL <- m2
exVars <- exVars.m2

# actually predict, on link scale so we can get proper CIs, exclude
treatPred <- as.data.frame(predict(MODEL, treatEff, type = 'link', se.fit = TRUE,
                                   exclude = exVars)) # <- pass in vars to exclude

## bind predictions to data we predicted at
treatPred <- cbind(treatEff, treatPred)

## extract inverse of link function; this is basically just `exp()`
ilink <- family(MODEL)$linkinv

## form 95% bayesian credible interval / frequentist across-function confidence interval
treatPred <- transform(treatPred, Fitted = ilink(fit),
                       Upper = ilink(fit + (2 * se.fit)),
                       Lower = ilink(fit - (2 * se.fit)))

## plot
p.plt <- ggplot(treatPred, aes(x = censusdate, y = Fitted)) +
    geom_point(data = rodent, mapping = aes(y = dipos, colour = treatment)) +
    geom_ribbon(aes(ymax = Upper, ymin = Lower, fill = treatment),
                alpha = 0.2) +
    geom_line(aes(colour = treatment)) +
    labs(y = 'Count', x = NULL) +
    theme(legend.position = 'top') +
    scale_colour_brewer(name = 'Treatment', type = 'qual', palette = 'Dark2') +
    scale_fill_brewer(name = 'Treatment', type = 'qual', palette = 'Dark2')
p.plt


## OK - figured out what the issue is. This is more the plot you should be looking at when you consider the *differences*, 
##    because the splines are fitted on the link scale and that's also the scale where differences make sense.
exVars2 <- c(exVars, "oTreatment")
treatPred2 <- predict(MODEL, treatEff, type = 'terms', exclude = exVars2)
treatPred2 <- rowSums(treatPred2)
treatPred2 <- cbind(treatEff, Fitted = treatPred2)
## plot
p.plt2 <- ggplot(treatPred2, aes(x = censusdate, y = Fitted)) +
    geom_line(aes(colour = treatment)) +
    labs(y = 'Count', x = NULL) +
    theme(legend.position = 'top') +
    scale_colour_brewer(name = 'Treatment', type = 'qual', palette = 'Dark2')
p.plt2
## Now compare the smooths on that plot with the differences I generate below

## Compute pairwise differences of smooths when fitted using ordered factors
osmooth_diff <- function(model, newdata, smooth_var, f1, f2, var, alpha = 0.05,
                         unconditional = FALSE, ...) {
    xp <- predict(model, newdata = newdata, type = 'lpmatrix', ...)
    ## reference level
    ref_level <- levels(newdata[[var]])[1L]
    ref_smooth <- grepl(paste0("s\\(", smooth_var, "\\)\\.{1}[[:digit:]]+$"), colnames(xp))
    c1 <- ref_smooth | grepl(f1, colnames(xp))
    c2 <- ref_smooth | grepl(f2, colnames(xp))
    r1 <- newdata[[var]] == f1
    r2 <- newdata[[var]] == f2
    ## difference rows of xp for data from comparison
    X <- xp[r1, ] - xp[r2, ]
    ## zero out cols of X related to splines for other lochs
    X[, ! (c1 | c2)] <- 0
    ## zero out the parametric cols
    X[, !grepl('^s\\(', colnames(xp))] <- 0
    dif <- X %*% coef(model)
    se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
    crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
    upr <- dif + (crit * se)
    lwr <- dif - (crit * se)
    data.frame(pair = paste(f1, f2, sep = '-'),
               diff = dif,
               se = se,
               upper = upr,
               lower = lwr)
}

## Compute pairwise diffs
d1 <- osmooth_diff(MODEL, treatEff, "numericdate", "CC", "EC", var = "oTreatment")
d2 <- osmooth_diff(MODEL, treatEff, "numericdate", "CC", "XC", var = "oTreatment")
d3 <- osmooth_diff(MODEL, treatEff, "numericdate", "EC", "XC", var = "oTreatment")
## stick together
diffs <- rbind(d1, d2, d3)
diffs <- with(treatEff, cbind(censusdate, diffs)) # bind on censusdate for plotting

## plot difference, facetted
ggplot(diffs, aes(x = censusdate, y = diff)) +
    geom_ribbon(aes(ymax = upper, ymin = lower), alpha = 0.2) +
    geom_line() +
    facet_wrap(~ pair, ncol = 1) +
    labs(y = 'Difference', x = NULL)

### Temporarily just use Erica's original data
dipoN_avg <- read_csv('Dipo_abundance_by_treatment.csv')
dipoN_avg <- mutate(dipoN_avg, oTreatment = ordered(treatment, levels = c('CC','EC','XC')))

## Add weights; to not change the log-likelihood we need normalized weights
dipoN_avg <- mutate(dipoN_avg,
                    wi = case_when(treatment %in% c('EC', 'XC') ~ 3,
                                   treatment == "CC" ~ 4),
                    wi = wi / mean(wi)) # normalize weights

## fit model with response conditionally ~ Tweedie, with weights = wi
## needed to increase k for all terms, but still some suggestion of residual pattern
m <- gam(dipos ~ oTreatment + s(numericdate, k = 25) + s(numericdate, by = oTreatment, k = 15),
         data = dipoN_avg, method='REML', family = tw, weights = wi)

gam.check(m, rep = 100)
summary(m)
plot(m, pages = 1, shade = TRUE, scale = 0)

## predict for treatment effects
pdataPooled <- with(dipoN_avg,
                    expand.grid(treatment = c('CC','EC','XC'),
                                censusdate = seq(min(censusdate), max(censusdate), length = np)))
pdataPooled <- mutate(pdataPooled,
                      oTreatment = ordered(treatment, levels = c('CC','EC','XC')),
                      numericdate = as.numeric(censusdate) / 1000)

## do predictions and created fitted and upper/lower limits of a CI
predPooled <- predict(m, newdata = pdataPooled, type = "link", se.fit = TRUE)
predPooled <- cbind(pdataPooled, as.data.frame(predPooled))
invTW <- family(m)$linkinv
predPooled <- mutate(predPooled, Fitted = invTW(fit),
                     Upper = invTW(fit + (2 * se.fit)),
                     Lower = invTW(fit - (2 * se.fit)))

## plot
tw.plt <- ggplot(predPooled, aes(x = censusdate, y = Fitted)) +
    geom_point(aes(y = dipos, x = censusdate, colour = treatment), data = dipoN_avg) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = treatment), alpha = 0.2) +
    geom_line(aes(colour = treatment)) +
    theme(legend.position = 'top') +
    scale_fill_brewer(name = 'Treatment', type = 'qual', palette = 'Dark2') +
    scale_colour_brewer(name = 'Treatment', type = 'qual', palette = 'Dark2') +
    labs(x = NULL, y = expression(Density ~ (ind. ~ plot^{-1})))
tw.plt

## compare the polt-level Poisson fit with the pooled Tweedie fit
comp.plt <- plot_grid(tw.plt + labs(title = "Pooled data; Tweedie response") + lims(y = c(0, 14)),
                      p.plt + labs(title = "Unpooled data; Poisson response") + lims(y = c(0, 14)),
                      nrow = 1, align = "hv")

comp.plt

ggsave("compare-models.pdf", comp.plt, height = 6, width = 12)
