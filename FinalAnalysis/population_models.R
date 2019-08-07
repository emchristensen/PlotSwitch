# PP Shifts Paper
# Ellen K. Bledsoe, with code from S. Supp
# March 7, 2018
# Modified by Erica Christensen for Plot Switch project 7/31/19
# This code runs population models for all species, though only D. merriami and D. ordii are discussed in the paper

# LIBRARIES and SOURCE CODE #

library(tidyverse)
library(RMark)


source("FinalAnalysis/population_model_functions.r")

# DATA FILES (PortalData folder was downloaded as part of create_data_series.R script)

# rodent file
rdat <- read.csv('PortalData/Rodents/Portal_rodent.csv', header = TRUE, na.strings = c(""), stringsAsFactors = FALSE)

# trapping data file
tdat <- read.csv('PortalData/Rodents/Portal_rodent_trapping.csv', header = T, na.strings = c(""), stringsAsFactors = F)

# plot treatments: type of switch
pdat <- data.frame(plot = 1:24, treatment = c('CC','CE','EE','CC','XC','EC',
                                              'XC','CE','CX','XX','CC','CX',
                                              'EC','CC','EE','XX','CC','EC',
                                              'EE','EE','EE','CE','XX','XC'))

##########################################################
# DATA PREP
##########################################################
# restrict to only the plots relevant to this project (controls after the switch in 2015)
plotswitchplots = c(4,11,14,17,6,13,18,5,7,24)

rdat_filtered = dplyr::filter(rdat, period>=437, plot %in% plotswitchplots)

#############################################################
# run RMARK models on each species of interest; save to csvs in Data/PopModelBest_afteronly
#############################################################
# this code runs population models and saves the output to csvs in Data/PopModelBest_afteronly
# the file names should be changed manually to reflect the date run

date_run = '20190731'
run_species_pop_model(rdat_filtered, sp='DM', date_run=date_run)
run_species_pop_model(rdat_filtered, sp='DO', date_run=date_run)
# small granivores 
run_species_pop_model(rdat_filtered, sp='RM', date_run=date_run)
run_species_pop_model(rdat_filtered, sp='PP', date_run=date_run)
run_species_pop_model(rdat_filtered, sp='BA', date_run=date_run)
run_species_pop_model(rdat_filtered, sp='PE', date_run=date_run)
# PF, PB, PL, PM, RF, RO have very few individuals
run_species_pop_model(rdat_filtered, sp='PF', date_run=date_run)
run_species_pop_model(rdat_filtered, sp='PB', date_run=date_run)
run_species_pop_model(rdat_filtered, sp='PM', date_run=date_run)
run_species_pop_model(rdat_filtered, sp='RO', date_run=date_run)


#############################################################
# make plots
#############################################################

# read in Mark results
dm_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_DM_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
do_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_DO_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
rm_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_RM_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
pp_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_PP_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
ba_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_BA_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
pe_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_PE_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
pm_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_PM_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
ro_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_RO_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
pf_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_PF_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
pb_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_PB_top_model_summary_",date_run,".csv"), stringsAsFactors = F)

# plot RMark results -- survival estimate
dm_plotdat <- prep_RMark_data_for_plotting(dm_results)
plot_estimated_survival(dplyr::filter(dm_plotdat,metric =='S'), paste0('D. merriami: n = ',dm_results$n_indiv[1]))

do_plotdat <- prep_RMark_data_for_plotting(do_results)
plot_estimated_survival(dplyr::filter(do_plotdat, metric=='S'), paste0('D. ordii: n = ',do_results$n_indiv[1]))

pb_plotdat <- prep_RMark_data_for_plotting(pb_results)
plot_estimated_survival(dplyr::filter(pb_plotdat, metric=='S'), paste0('C. baileyi: n = ',pb_results$n_indiv[1]))

pp_plotdat <- prep_RMark_data_for_plotting(pp_results)
plot_estimated_survival(dplyr::filter(pp_plotdat, metric=='S'), paste0('C. penicillatus: n = ',pp_results$n_indiv[1]))

ba_plotdat <- prep_RMark_data_for_plotting(ba_results)
plot_estimated_survival(dplyr::filter(pp_plotdat, metric=='S'), paste0('B. taylori: n = ',ba_results$n_indiv[1]))

pf_plotdat <- prep_RMark_data_for_plotting(pf_results)
plot_estimated_survival(dplyr::filter(pf_plotdat, metric=='S'), paste0('P. flavus: n = ',pf_results$n_indiv[1]))

pe_plotdat <- prep_RMark_data_for_plotting(pe_results)
plot_estimated_survival(dplyr::filter(pe_plotdat, metric=='S'), paste0('P. eremicus: n = ',pe_results$n_indiv[1]))

pm_plotdat <- prep_RMark_data_for_plotting(pm_results)
plot_estimated_survival(dplyr::filter(pm_plotdat, metric=='S'), paste0('P. maniculatus: n = ',pm_results$n_indiv[1]))

ro_plotdat <- prep_RMark_data_for_plotting(ro_results)
plot_estimated_survival(dplyr::filter(ro_plotdat, metric=='S'), paste0('R. montanus: n = ',ro_results$n_indiv[1]))

rm_plotdat <- prep_RMark_data_for_plotting(rm_results)
plot_estimated_survival(dplyr::filter(rm_plotdat, metric=='S'), paste0('R. megalotis: n = ',rm_results$n_indiv[1]))


# plot Rmark results -- transition (psi) only DM and RM have difference in psi by treatment, and RM doesn't have enough info
plot_estimated_survival(dplyr::filter(dm_plotdat, metric=='Psi'), paste0('D. merriami: n = ',dm_results$n_indiv[1]))
#plot_estimated_survival(dplyr::filter(rm_plotdat, metric=='Psi'), paste0('R. megalotis: n = ',rm_results$n_indiv[1]))

#############################################################
# Number of New Individuals Showing Up on Plots
#############################################################
# run this code for both "DM" and "DO"
new_per_plot_trt = new_captures_by_plot('DO',rdat,tdat)

# remove rows where plots were not sampled
new_per_plot_trt = new_per_plot_trt[new_per_plot_trt$sampled==1,]
new_per_plot_trt$yr = new_per_plot_trt$year
new_per_plot_trt$yr[new_per_plot_trt$month<4] <- new_per_plot_trt$yr[new_per_plot_trt$month<4]-1

# create column for "years after switch"
new_per_plot_trt$yr_since_change <- NA
new_per_plot_trt$yr_since_change[new_per_plot_trt$period %in% 437:447] <- 1
new_per_plot_trt$yr_since_change[new_per_plot_trt$period %in% 448:460] <- 2
new_per_plot_trt$yr_since_change[new_per_plot_trt$period %in% 461:472] <- 3
new_per_plot_trt$yr_since_change[new_per_plot_trt$period > 472] <- 4 # April-July 2018

# calculate total new animals per treatment type in each year since switch
t =new_per_plot_trt %>% group_by(yr_since_change, treatment) %>%
  summarise(total_new_per_year=sum(count)) %>%
  spread(yr_since_change, total_new_per_year)
# WARNING: there are 4 CC plots and only 3 of each EC and XC plots
# number of new animals per plot per year: (these are the numbers reported in the supplement)
cbind(t[,1],t[,2:4]/c(4,3,3))


##########################################################################################
# scraps
##############################################################################################

# ==============================================================
# another way of looking at new individuals arriving on plots
# ==============================================================
# 
# new_per_plot_trt = new_captures_by_plot('DM',rdat,tdat)
# 
# # remove rows where plots were not sampled
# new_per_plot_trt = new_per_plot_trt[new_per_plot_trt$sampled==1,]
# new_per_plot_trt$yr = new_per_plot_trt$year
# new_per_plot_trt$yr[new_per_plot_trt$month<4] <- new_per_plot_trt$yr[new_per_plot_trt$month<4]-1
# 
# 
# new_per_plot_yr = aggregate(new_per_plot_trt$count, by=list(plot=new_per_plot_trt$plot,
#                                                             treatment=new_per_plot_trt$treatment,
#                                                             yr=new_per_plot_trt$yr),
#                             FUN=sum) %>% rename(count=x) %>% filter(yr<2018)
# 
# test = data.frame(treatment=new_per_plot_yr$treatment, count=new_per_plot_yr$count, yr=new_per_plot_yr$yr)
# 
# t = new_per_plot_yr %>% group_by(treatment, yr) %>%
#   summarise(total=sum(count))
# 
# # plot: total individuals by treatment and year. boxplots represent spread from multiple plots (3-4)
# ggplot(test, aes(x=treatment, y=count, fill=treatment)) +
#   geom_boxplot() +
#   geom_point(position=position_jitter(.2)) +
#   ylab('# new individuals per plot per year') +
#   ggtitle('DO')
# 
# 
# # code for running MARK model including before/after switch data, with a factor for before/after
# rdat_filtered = dplyr::filter(rdat, period>=415, plot %in% c(4,11,14,17,6,13,18,5,7,24))
# dmdat = individual_tag_cleanup(c('DM','DO'), rdat_filtered)
# dmdat_trt = merge(dmdat, pdat, by=c('plot'))
# 
# mark_dm = create_trmt_hist(dmdat_trt)
# # resulting warnings are probably for animals captured twice in same sampling period
# 
# # prep data for RMark
# all_ms <- select(mark_dm, captures) %>% dplyr::rename("ch" = "captures")
# first_per <- min(dmdat_trt$period)
# 
# # Process data
# ms.pr = process.data(all_ms, begin.time = first_per, model = "Multistrata")
# 
# # Create default design data
# ms.ddl = make.design.data(ms.pr)
# 
# # add design covariates for PB era
# after_switch = as.factor(437:476)
# 
# ms.ddl$S$after_switch = 0
# ms.ddl$S$after_switch[ms.ddl$S$time %in% after_switch] = 1
# 
# ms.ddl$p$after_switch = 0
# ms.ddl$p$after_switch[ms.ddl$p$time %in% after_switch] = 1
# 
# ms.ddl$Psi$after_switch = 0
# ms.ddl$Psi$after_switch[ms.ddl$Psi$time %in% after_switch] = 1
# 
# # Run the models and examine the output
# ms.results = run.ms(S_dot = NULL,
#                     S_stratum = list(formula = ~ -1 + stratum * after_switch),
#                     p_dot = list(formula = ~ 1),
#                     p_stratum = NULL,
#                     Psi_s = list(formula =  ~ -1 + stratum:tostratum * after_switch, link = "logit"))
# rmark_results <- ms.results$S.stratum.p.dot.Psi.s$results$real
# rmark_results
# 
# # prep for plotting
# plot_rmark <- prep_RMark_data_for_plotting(rmark_results)
