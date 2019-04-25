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

# species file
#sdat <- read.csv('PortalData/Rodents/Portal_rodent_species.csv', header = TRUE, na.strings = c(""), stringsAsFactors = FALSE)

# trapping data
#tdat <- read.csv('PortalData/Rodents/Portal_rodent_trapping.csv', header = TRUE, stringsAsFactors = FALSE)

# plot treatments: type of switch
#pdat <- read.csv('PortalData/SiteandMethods/Portal_plots.csv', header = TRUE, stringsAsFactors = FALSE)
pdat <- data.frame(plot = 1:24, treatment = c('CC','CE','EE','CC','XC','EC',
                                              'XC','CE','CX','XX','CC','CX',
                                              'EC','CC','EE','XX','CC','EC',
                                              'EE','EE','EE','CE','XX','XC'))

##########################################################
# DATA PREP
##########################################################

#---------------------------------------------------------
# Clean the Data
#---------------------------------------------------------
# restrict to only the plots relevant to this project (controls after the switch in 2015)
rdat_filtered = dplyr::filter(rdat, period>=437, plot %in% c(4,11,14,17,6,13,18,5,7,24))

#dipodat = individual_tag_cleanup(sp=c('DO','DM','DS'), rdat_filtered)
#smgrandat = individual_tag_cleanup(sp=c('BA','PB','PP','PF','PE','PL','PM','RF','RM','RO'), rdat_filtered)


#############################################################
# DM models
#############################################################
dmdat = individual_tag_cleanup(c('DM'), rdat_filtered)
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
 
# run model just on after-switch data
ms.results = run.ms(S_dot = NULL,
                    S_stratum = list(formula = ~ -1 + stratum),
                    p_dot = list(formula = ~ 1),
                    p_stratum = NULL,
                    Psi_s = list(formula =  ~ -1 + stratum:tostratum, link = "logit"))
ms.results

rmark_results <- ms.summary$results$real
rmark_results
write.csv(ms.summary$results$real, "Data/MARK_DM_top_model_summary_[DATE]2.csv")

#############################################################
# PB models
#############################################################
pbdat = individual_tag_cleanup('PB', rdat_filtered)
pbdat_trt = merge(pbdat, pdat, by=c('plot'))

mark_pb = create_trmt_hist(pbdat_trt)
# resulting warnings are probably for animals captured twice in same sampling period

# prep data for RMark
all_ms <- select(mark_pb, captures) %>% dplyr::rename("ch" = "captures")
first_per <- min(pbdat_trt$period)

# Process data
ms.pr = process.data(all_ms, begin.time = first_per, model = "Multistrata")

# Create default design data
ms.ddl = make.design.data(ms.pr)

# add design covariates for before/after switch
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
ms.results
names(ms.results)

ms.summary = ms.results$S.stratum.p.dot.Psi.s

rmark_results <- ms.summary$results$real
rmark_results
write.csv(ms.summary$results$real, "Data/MARK_PB_top_model_summary_[DATE].csv")

##########################
# make plots
#######################

# read in Mark results if skipping that section
rmark_results <- read.csv("../PP_shifts/data/top_model_summary_20190416.csv", stringsAsFactors = FALSE)

# plot RMark results

plot_rmark <- prep_RMark_data_for_plotting(rmark_results)

(plot2a <- plot_estimated_survival(plot_rmark))

(plot2b <- plot_transition_probability(plot_rmark))

#------------------------------------------------------------
# Number of New PP Individuals Showing Up on Plots
#------------------------------------------------------------

# make empty dataframe
first_period <- setNames(data.frame(matrix(ncol = 21, nrow = 0)), names(PP_only))

# create dataframe of only first period each tag is present
for (i in 1:length(tags_all)){
  tmp <- PP_only[PP_only$tag == tags_all[i],] # rows for a given tag
  tmp2 <- tmp[tmp$period == min(tmp$period),] # row with earliest period
  first_period <- rbind(first_period, tmp2)   # add to new dataframe
}

# total number new PPs (avg plot sum by year)
new_PP_per_plot <- first_period %>%
  filter(plot_type != "Removal") %>%
  group_by(plot, month, year, plot_type) %>%
  summarise(count = n(species)) %>%
  ungroup()

new_PP_per_plot <- new_PP_per_plot %>% 
  group_by(year, plot_type) %>% 
  summarise(sum_by_year = sum(count))

new_PP_per_plot$time_point <- NA

for (i in 1:nrow(new_PP_per_plot)) {
  if (new_PP_per_plot$year[i] <= 1997) {
    new_PP_per_plot$time_point[i] = "Before"
  } else {
    new_PP_per_plot$time_point[i] = "After"
  }
}

anova2w <- aov(sum_by_year ~ plot_type * time_point, data = new_PP_per_plot)
summary(anova2w)

# plot new PP individuals
(plot2c <- plot_new_PP_individuals(new_PP_per_plot))

# Make Figure 2
(plot2 <- plot2a + plot2b - plot2c + plot_layout(ncol = 1))

# ggsave("figures/1989-2010/Figure2.png", plot2, width = 6, height = 7, dpi = 600)
# ggsave("figures/1989-2010/Figure2.tiff", plot2,
#        width = 6, height = 7, dpi = 600, compression = "lzw")


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

# Run the models and examine the output
# ms.results = run.ms(S_dot = NULL,
#                     S_stratum = list(formula = ~ -1 + stratum * after_switch),
#                     p_dot = list(formula = ~ 1),
#                     p_stratum = NULL,
#                     Psi_s = list(formula =  ~ -1 + stratum:tostratum * after_switch, link = "logit"))
# rmark_results <- ms.summary$results$real
# rmark_results