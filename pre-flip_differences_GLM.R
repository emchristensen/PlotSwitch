# comparing values (number of dipos, number of non-dipo species, and total energy) between treatments, before treatment change in 2015

library(multcomp)
library(dplyr)
library('mgcv')
library(ggplot2)

# dipos -----
df1 = read.csv('Dipo_counts.csv') %>% filter(numericdate<16.543)
# plot data
ggplot(df1,aes(x=censusdate,y=n,colour=treatment)) + geom_jitter(width=.1,height=.2)

# glm of krat abundance
dipo.glm= glm(n ~ treatment + numericdate, data = df1, family = poisson())
summary(dipo.glm)
confint(dipo.glm)

# test for pairwise differences using glht (general linear hypothesis)
summary(glht(dipo.glm, linfct=mcp(treatment='Tukey')))


# small granivores ----
df2 = read.csv('SmallGranivores.csv') %>% filter(numericdate<16.543)
ggplot(df2,aes(x=censusdate,y=n,colour=treatment)) + geom_jitter(width=.1,height=.2)

sm.glm = glm(n ~ treatment+numericdate, data = df2, family = poisson())
summary(sm.glm)
confint(sm.glm)
layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(sm.glm)
layout(1)

summary(glht(sm.glm, linfct=mcp(treatment='Tukey')))

# energy ---- 
df3 = read.csv('TotalCommunityEnergy.csv') %>% filter(numericdate<16.543)
ggplot(df3,aes(x=censusdate,y=n,colour=treatment)) +geom_jitter(width=.1,height=.2)

en.mod = gam(n ~ treatment+numericdate, data = df3, family = tw, method = "REML")
summary(en.mod)     
## confint(en.mod)                         # nope, fails for a gam

## mgcv's model.matrix.gam doesn't add the attributes that we
## need for multcom::glth(), but it does return the required
## objects in the fitted model:
##   This adds `contrasts` and `assign` attributes to the
##   model matrix returned by model.matrix.gam
`model.matrix.gam` <- function (object, ...) {
    if (!inherits(object, "gam")) 
        stop("`object' is not of class \"gam\"")
    out <- predict(object, type = "lpmatrix", ...)
    attr(out, "contrasts") <- object[["contrasts"]]
    attr(out, "assign") <- object[["assign"]]
    out
}

## as this is all in the mgxcv namespace we need to source the above
## function into R and then use the following to assign the value of
## the function into the mcgv namespace
assignInNamespace('model.matrix.gam', model.matrix.gam, ns = 'mgcv')

## estimate the Tukey all pairwise comparison of treatment levels
en.glht <- glht(en.mod, linfct=mcp(treatment='Tukey'))
summary(en.glht)                        # summary

## check I didn't mess up
newd <- expand.grid(treatment = c('CC', 'EC', 'XC'),
                    numericdate = 15.777)
pred <- predict(en.mod, newd)
pred[c(2,3,3)] - pred[c(1,1,2)] # should be Estimate column in summary(en.glht)
