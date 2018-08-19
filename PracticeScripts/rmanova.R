library(nlme)
library(car)
library(rcompanion)
library(dplyr)
library(lsmeans)
library(ggpubr)

source('data_functions.R')



#' @description prepare data for doing rmanova
#' @param df data frame containing censusdate and plot as columns
#' 
prepare_dataframe = function(df) {
  df$censusdate = as.Date(df$censusdate)
  timecol = data.frame(censusdate = unique(df$censusdate))
  timecol$ordereddate = as.numeric(row.names(timecol))
  dat = merge(df,timecol)
  df$plot = as.factor(df$plot)
  return(dat)
}


#===========================================================
# Just an anova on before-switch data, to show if there are differences in treatment before any changes
df = read.csv('Dipo_counts.csv')
df = read.csv('SpeciesRichness.csv')
df = read.csv('SmallGranivores.csv')
df = read.csv('TotalCommunityEnergy.csv')

dat = prepare_dataframe(df) %>% filter(ordereddate<22)

model = lme(n ~ treatment + ordereddate + treatment*ordereddate,
            random=~1|plot,
            data=dat,
            method='REML')

# post hoc
marginal = lsmeans(model, 
                   ~ treatment:ordereddate)

cld(marginal,
    alpha   = 0.05, 
    Letters = letters,     ### Use lower-case letters for .group
    adjust  = "tukey")


#' Attempting repeated measures anova
#' 
#' I want to compare the effect of putting krats back on removal plots
#' 
#' used http://rcompanion.org/handbook/I_09.html
#' 
#' 
#' Metrics:
#'   - number of Dipodomys spp
#'   - number of small granivores
#'   - number of Onychomys (unrelated group)
#' 
#' Choices:
#'   - 3 month avg for figures only
#'   - corAR1 (time variable must be an integer variable)
#'

# data viz
# dataSum = rcompanion::groupwiseMean(n ~ treatment + numericdate,
#                     data   = smgran,
#                     conf   = 0.95,
#                     digits = 3,
#                     traditional = FALSE,
#                     percentile  = TRUE)
# 
# dataSum
# 
# pd = position_dodge(.2)
# 
# ggplot2::ggplot(dataSum, aes(x =    numericdate,
#                 y =    Mean,
#                 color = treatment)) +
#   geom_errorbar(aes(ymin=Percentile.lower,
#                     ymax=Percentile.upper),
#                 width=.2, size=0.7, position=pd) +
#   geom_point(shape=15, size=4, position=pd) +
#   theme_bw() +
#   theme(axis.title = element_text(face = "bold")) +
#   ylab("Mean small granivores")
# 



# attempting rmanova
smgran = read.csv('SmallGranivores.csv')
smgran$censusdate = as.Date(smgran$censusdate)

#quarterly = avg_3month(smgran)
# make an integer time column
timecol = data.frame(censusdate =unique(smgran$censusdate))
timecol$ordereddate = as.numeric(row.names(timecol))

smgran = merge(smgran,timecol)
smgran$plot = as.factor(smgran$plot)

smgranbefore = filter(smgran,ordereddate<22)

# figure out autocorrelation structure -- they took the lag1 values from this and put it in the next models
# model with no random factor
model.a = gls(n ~ treatment + ordereddate + treatment*ordereddate,
              data=smgranbefore)
ACF(model.a,
    form = ~ ordereddate | plot)

# model with plot as random factor
model.b = lme(n ~ treatment + ordereddate + treatment*ordereddate,
              random=~1|plot,
              data=smgranbefore)
ACF(model.b)

# ======================================================
# take a stab at a model
model = lme(n ~ treatment + ordereddate + treatment*ordereddate,
            random=~1|plot,
           # correlation = corAR1(form=~ordereddate | plot, value = 0.349),
            data=smgranbefore,
            method='REML')
Anova(model)

model.fixed = gls(n ~ treatment + ordereddate + treatment*ordereddate,
            #      corAR1(form = ~ordereddate|plot, value = 0.349),
            data=smgranbefore,
            method='REML')

anova(model,model.fixed)
# I think this is saying that the model doesn't need the random effect of plot

# compare to a null model
model.null = lme(n ~ 1,
                 random = ~1|plot,
                 data = smgranbefore)

nagelkerke(model, 
           model.null)

# post hoc test
library(lsmeans)

marginal = lsmeans(model, 
                   ~ treatment:ordereddate)

cld(marginal,
    alpha   = 0.05, 
    Letters = letters,     ### Use lower-case letters for .group
    adjust  = "tukey")


