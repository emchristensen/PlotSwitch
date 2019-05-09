# PP Shifts Paper
# Ellen K. Bledsoe, with code from S. Supp
# March 7, 2018
# Modified EMC for Plot Switch project 4/22/19

# LIBRARIES and SOURCE CODE #

library(tidyverse)
#library(portalr)
#library(plotrix)
library(RMark)
#library(forecast)
#library(nlme)
#library(patchwork) # devtools::install_github("thomasp85/patchwork")
#library(rapportools)

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

#rdat_filtered = dplyr::filter(rdat, period>=437, plot %in% plotswitchplots)
rdat_filtered = dplyr::filter(rdat, period>=413, plot %in% plotswitchplots)


#############################################################
# run RMARK models on each species of interest; save to csv
#############################################################
# this code runs population models and saves the output to csvs
# the file names should be changed manually to reflect the date run
run_species_pop_model(rdat_filtered, sp='DM')
run_species_pop_model(rdat_filtered, sp='DO')
# small granivores
run_species_pop_model(rdat_filtered, sp='PB')
run_species_pop_model(rdat_filtered, sp='PP')
run_species_pop_model(rdat_filtered, sp='BA')
run_species_pop_model(rdat_filtered, sp='PF')
run_species_pop_model(rdat_filtered, sp='PE')
run_species_pop_model(rdat_filtered, sp='PL')
run_species_pop_model(rdat_filtered, sp='PM')
run_species_pop_model(rdat_filtered, sp='RF')
run_species_pop_model(rdat_filtered, sp='RO')
run_species_pop_model(rdat_filtered, sp='RM')

#############################################################
# make plots
#############################################################

# read in Mark results
dm_results <- read.csv("Data/MARK_DM_top_model_summary_[DATE].csv", stringsAsFactors = F)
do_results <- read.csv("Data/MARK_DO_top_model_summary_[DATE].csv", stringsAsFactors = F)
pb_results <- read.csv("Data/MARK_PB_top_model_summary_[DATE].csv", stringsAsFactors = F)
pp_results <- read.csv("Data/MARK_PP_top_model_summary_[DATE].csv", stringsAsFactors = F)
ba_results <- read.csv("Data/MARK_BA_top_model_summary_[DATE].csv", stringsAsFactors = F)
pf_results <- read.csv("Data/MARK_PF_top_model_summary_[DATE].csv", stringsAsFactors = F)
pe_results <- read.csv("Data/MARK_PE_top_model_summary_[DATE].csv", stringsAsFactors = F)
pl_results <- read.csv("Data/MARK_PL_top_model_summary_[DATE].csv", stringsAsFactors = F)
pm_results <- read.csv("Data/MARK_PM_top_model_summary_[DATE].csv", stringsAsFactors = F)
rf_results <- read.csv("Data/MARK_RF_top_model_summary_[DATE].csv", stringsAsFactors = F)
ro_results <- read.csv("Data/MARK_RO_top_model_summary_[DATE].csv", stringsAsFactors = F)
rm_results <- read.csv("Data/MARK_RM_top_model_summary_[DATE].csv", stringsAsFactors = F)

# plot RMark results
dm_plotdat <- prep_RMark_data_for_plotting(dm_results)
plot_estimated_survival(dm_plotdat, paste0('D. merriami: n = ',dm_results$n_indiv[1]))

do_plotdat <- prep_RMark_data_for_plotting(do_results)
plot_estimated_survival(do_plotdat, paste0('D. ordii: n = ',do_results$n_indiv[1]))

pb_plotdat <- prep_RMark_data_for_plotting(pb_results)
plot_estimated_survival(pb_plotdat, paste0('C. baileyi: n = ',pb_results$n_indiv[1]))

pp_plotdat <- prep_RMark_data_for_plotting(pp_results)
plot_estimated_survival(pp_plotdat, paste0('C. penicillatus: n = ',pp_results$n_indiv[1]))

ba_plotdat <- prep_RMark_data_for_plotting(ba_results)
plot_estimated_survival(ba_plotdat, paste0('B. taylori: n = ',ba_results$n_indiv[1]))

pf_plotdat <- prep_RMark_data_for_plotting(pf_results)
plot_estimated_survival(pf_plotdat, paste0('P. flavus: n = ',pf_results$n_indiv[1]))

pe_plotdat <- prep_RMark_data_for_plotting(pe_results)
plot_estimated_survival(pe_plotdat, paste0('P. eremicus: n = ',pe_results$n_indiv[1]))

pl_plotdat <- prep_RMark_data_for_plotting(pl_results)
plot_estimated_survival(pl_plotdat, paste0('P. leucopus: n = ',pl_results$n_indiv[1]))

pm_plotdat <- prep_RMark_data_for_plotting(pm_results)
plot_estimated_survival(pm_plotdat, paste0('P. maniculatus: n = ',pm_results$n_indiv[1]))

rf_plotdat <- prep_RMark_data_for_plotting(rf_results)
plot_estimated_survival(rf_plotdat, paste0('R. fulvescens: n = ',rf_results$n_indiv[1]))

ro_plotdat <- prep_RMark_data_for_plotting(ro_results)
plot_estimated_survival(ro_plotdat, paste0('R. montanus: n = ',ro_results$n_indiv[1]))

rm_plotdat <- prep_RMark_data_for_plotting(rm_results)
plot_estimated_survival(rm_plotdat, paste0('R. megalotis: n = ',rm_results$n_indiv[1]))

#############################################################
# Number of New Individuals Showing Up on Plots
#############################################################
new_per_plot_trt = new_captures_by_plot('DM',rdat,tdat)

# remove rows where count = NA: these are plots that were not sampled
new_per_plot_trt = new_per_plot_trt[!is.na(new_per_plot_trt$count),]
# get average by treatment

new_per_trt = aggregate(new_per_plot_trt$count, by=list(period=new_per_plot_trt$period, treatment=new_per_plot_trt$treatment),
                        FUN=mean)

summarize(group_by(new_per_trt, treatment),
          mean=mean(x), sd=sd(x), min=min(x), max=max(x))

ggplot(new_per_trt, aes(x=period,y=x,colour=treatment)) +
  geom_point() +
  geom_line()

anova2w <- aov(x ~ treatment * period, data = new_per_trt)
summary(anova2w)
anova2w <- aov(x ~ treatment, data=new_per_trt)
summary(anova2w)

# plot new PP individuals
#(plot2c <- plot_new_PP_individuals(new_PP_per_plot))

# Make Figure 2
#(plot2 <- plot2a + plot2b - plot2c + plot_layout(ncol = 1))


### RUN ANOVAS ###
# Things to deal with:
#   - which way to summarize the data
#   - what about after 2010? seems to be washing out treatment and interaction

# 2-way ANOVA by monthly count

# new_PP_per_plot$time_point <- NA
# 
# for (i in 1:nrow(new_PP_per_plot)) {
#   if (new_PP_per_plot$year[i] < 1997) {
#     new_PP_per_plot$time_point[i] = "Before"
#   } else if (new_PP_per_plot$year[i] > 1997){
#     new_PP_per_plot$time_point[i] = "After"
#   } else {
#     if (new_PP_per_plot$month[i] < 7){
#       new_PP_per_plot$time_point[i] = "Before"
#     } else {
#       new_PP_per_plot$time_point[i] = "After"
#     }
#   }
# }
# 
# new_PP_per_plot <- filter(new_PP_per_plot, year <= 2010)
# anova2w <- aov(count ~ plot_type * time_point, data = new_PP_per_plot)
# summary(anova2w)
# 


# 2-way ANOVA by year

# new_PP_per_plot <- new_PP_per_plot %>% 
#   group_by(year, plot_type) %>% 
#   summarise(count = sum(count))
# 
# new_PP_per_plot$time_point <- NA
# 
# for (i in 1:nrow(new_PP_per_plot)) {
#   if (new_PP_per_plot$year[i] <= 1997) {
#     new_PP_per_plot$time_point[i] = "Before"
#   } else {
#     new_PP_per_plot$time_point[i] = "After"
#   }
# }
# 
# new_PP_per_plot <- filter(new_PP_per_plot, year <= 2010)
# anova2w <- aov(count ~ plot_type * time_point, data = new_PP_per_plot)
# summary(anova2w)
# 
# ggplot(data = new_PP_per_plot, aes(x = year, y = count, color = plot_type)) +
#   geom_point() +
#   geom_line() +
#   theme_bw()

# 2-way ANOVA by avg per plot per year

# new_PP_per_plot_summary$time_point <- NA
# 
# for (i in 1:nrow(new_PP_per_plot_summary)) {
#   if (new_PP_per_plot_summary$year[i] <= 1997) {
#     new_PP_per_plot_summary$time_point[i] = "Before"
#   } else {
#     new_PP_per_plot_summary$time_point[i] = "After"
#   }
# }
# 
# new_PP_per_plot_summary <- filter(new_PP_per_plot_summary, year <= 2010)
# anova2w <- aov(avg_plot_sum_by_year ~ plot_type * time_point, data = new_PP_per_plot_summary)
# summary(anova2w)

#############################################################
# System-level Aspects of Patch Preference
#############################################################

# download biomass data by plot from portalr
biomass_data <- portalr::biomass(path = "repo", level = "Plot")
energy_data <- portalr::energy(path = "repo", level = "Plot")

# select certain treatments and filter by time
biomass_dat <- biomass_data %>%
  filter(treatment == "control" | treatment == "exclosure", 
         period >= 118 & period <= 433) # get the right time periods
energy_dat <- energy_data %>%
  filter(treatment == "control" | treatment == "exclosure", 
         period >= 118 & period <= 433) 


# add a year column for later summarization
year_prd_pairs <- unique(tdat[,c("year", "period")]) # get associated years and periods
biomass_dat$year = NA
energy_dat$year = NA

for (i in 1:nrow(biomass_dat)){
  prd <- biomass_dat$period[i]
  biomass_dat$year[i] = year_prd_pairs$year[year_prd_pairs$period == prd]
}

for (i in 1:nrow(energy_dat)){
  prd <- energy_dat$period[i]
  energy_dat$year[i] = year_prd_pairs$year[year_prd_pairs$period == prd]
}

#------------------------------------------------------------
# Biomass Ratios
#------------------------------------------------------------

# sum across rows and rename column
biomass_dat_rowSums <- as.data.frame(rowSums(biomass_dat[,4:24]))
colnames(biomass_dat_rowSums) <- c("rowSums")

# summarise biomass to get total by period and plot type
biomass_total <- cbind(biomass_dat, biomass_dat_rowSums) %>%
  group_by(year, treatment) %>%
  summarise(totals = sum(rowSums))

# change the data structure to run the linear model
biomass_spread <- tidyr::spread(biomass_total, treatment, totals)

# ratio
biomass_ratio <- biomass_spread %>% mutate(EX_to_CO_ratio = exclosure/control)

(plot3 <- plot_biomass_ratio(biomass_ratio))
# ggsave("figures/1989-2010/Figure3.png", plot3, width = 3.5, height = 3, dpi = 600)

## ENERGY ##

# sum across rows and rename column
energy_dat_rowSums <- as.data.frame(rowSums(energy_dat[,4:24]))
colnames(energy_dat_rowSums) <- c("rowSums")

# summarise energy to get total by period and plot type
energy_total <- cbind(energy_dat, energy_dat_rowSums) %>%
  group_by(year, treatment) %>%
  summarise(totals = sum(rowSums))

# change the data structure to run the linear model
energy_spread <- tidyr::spread(energy_total, treatment, totals)

# ratio
energy_ratio <- energy_spread %>% mutate(EX_to_CO_ratio = exclosure/control)

(plot3_energy <- plot_energy_ratio(energy_ratio))
# ggsave("figures/1989-2010/Figure3_energy.png", plot3_energy, width = 3.5, height = 3, dpi = 600)

##########################################################################################
# scraps
##############################################################################################

# code for running MARK model including before/after switch data, with a factor for before/after
rdat_filtered = dplyr::filter(rdat, period>=415, plot %in% c(4,11,14,17,6,13,18,5,7,24))
dmdat = individual_tag_cleanup(c('DM','DO'), rdat_filtered)
dmdat_trt = merge(dmdat, pdat, by=c('plot'))

mark_dm = create_trmt_hist(dmdat_trt)
# resulting warnings are probably for animals captured twice in same sampling period

# prep data for RMark
all_ms <- select(mark_dm, captures) %>% dplyr::rename("ch" = "captures")
first_per <- min(dmdat_trt$period)

# Process data
ms.pr = process.data(all_ms, begin.time = first_per, model = "Multistrata")

# Create default design data
ms.ddl = make.design.data(ms.pr)

# add design covariates for PB era
after_switch = as.factor(437:476)

ms.ddl$S$after_switch = 0
ms.ddl$S$after_switch[ms.ddl$S$time %in% after_switch] = 1

ms.ddl$p$after_switch = 0
ms.ddl$p$after_switch[ms.ddl$p$time %in% after_switch] = 1

ms.ddl$Psi$after_switch = 0
ms.ddl$Psi$after_switch[ms.ddl$Psi$time %in% after_switch] = 1

# Run the models and examine the output
ms.results = run.ms(S_dot = NULL,
                    S_stratum = list(formula = ~ -1 + stratum * after_switch),
                    p_dot = list(formula = ~ 1),
                    p_stratum = NULL,
                    Psi_s = list(formula =  ~ -1 + stratum:tostratum * after_switch, link = "logit"))
rmark_results <- ms.results$S.stratum.p.dot.Psi.s$results$real
rmark_results

# prep for plotting
plot_rmark <- prep_RMark_data_for_plotting(rmark_results)
