# This code runs population models for D. merriami and D. ordii
# EMC 7/2019

library(tidyverse)
library(RMark)
library(cowplot)

source("FinalAnalysis/population_model_functions.r")
theme_set(theme_bw())
cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3", "#0072B2", "#D55E00", 
                "#CC79A7", "#F0E442")

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

#############################################################
# make plots
#############################################################
# read in Mark results
dm_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_DM_top_model_summary_",date_run,".csv"), stringsAsFactors = F)
do_results <- read.csv(paste0("Data/PopModelBest_afteronly/MARK_DO_top_model_summary_",date_run,".csv"), stringsAsFactors = F)

# plot RMark results -- survival estimate
dm_plotdat <- prep_RMark_data_for_plotting(dm_results)
dm_survival_plot <- plot_estimated_survival(dplyr::filter(dm_plotdat,metric =='S'), colorvalues=cbbPalette[c(6,1,4)], maintitle='')

do_plotdat <- prep_RMark_data_for_plotting(do_results)
do_survival_plot <- plot_estimated_survival(dplyr::filter(do_plotdat, metric=='S'), colorvalues=cbbPalette[c(6,1,4)], maintitle='')

# plot Rmark results -- transition (psi) only DM has significant difference in psi by treatment
dm_transitiondat = dplyr::filter(dm_plotdat, metric=='Psi') %>% 
  mutate(Transition=Treatment,
         Treatment=c('Control','Control','Kangaroo rat+','Kangaroo rat+','Rodent+','Rodent+'))
dm_transition_plot <- plot_estimated_survival(dm_transitiondat, colorvalues=cbbPalette[c(6,1,4)], maintitle='') + 
  theme(legend.position = 'top',
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.margin = margin(t=0, unit='cm'))

#############################################################
# Number of New Individuals Showing Up on Plots
#############################################################
# run function to get number of new individuals of given species per year
new_per_year_trt_do <- new_captures_by_year('DO',rdat,tdat)

# group by year, color by treatment: DO
ggplot(new_per_year_trt_do, aes(x=xposition, y=total_new_per_year, color=treatment)) +
  geom_jitter(alpha=.5, width=.03) +
  geom_errorbarh(aes(xmin=xposition-.05, xmax=xposition+.05, y=mean_by_trt, color=treatment), size=3.5) +
  ylab('Average # new per year') +
  xlab('') +
  scale_x_continuous(labels=c('Year 1','Year 2','Year 3'), breaks=c(1.2,2.2,3.2)) +
  scale_color_manual(name='Treatment:',
                     breaks=c('CC','EC','XC'),
                     labels=c('Control','Kangaroo rat+','Rodent+'),
                     values = cbbPalette[c(6,1,4)]) +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        plot.margin = margin(l = 5, t = 0, b = 0))

new_per_year_trt_dm <- new_captures_by_year('DM', rdat, tdat)

# group by year, color by treatment: DM
dm_newplot <- ggplot(new_per_year_trt_dm, aes(x=xposition, y=total_new_per_year, color=treatment)) +
  geom_jitter(alpha=.5, width=.03) +
  geom_errorbarh(aes(xmin=xposition-.05, xmax=xposition+.05, y=mean_by_trt, color=treatment), size=3.5) +
  ylab('Average # new per year') +
  xlab('') +
  scale_x_continuous(labels=c('Year 1','Year 2','Year 3'), breaks=c(1.2,2.2,3.2)) +
  scale_color_manual(name='Treatment:',
                    breaks=c('CC','EC','XC'),
                    labels=c('Control','Kangaroo rat+','Rodent+'),
                    values = cbbPalette[c(6,1,4)]) +
  theme(axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 6),
        legend.position = "none",
        plot.margin = margin(l = 5, t = 0, b = 0))
dm_newplot
##########################################################
# Figures for manuscript
##########################################################

## Cowplot grid
dm_pop_comboplot = plot_grid(dm_newplot, dm_survival_plot, labels = c('','B'), ncol = 2, align='h')
dm_pop_comboplot
dm_pop_comboplot2 = plot_grid(dm_pop_comboplot, NULL, dm_transition_plot, rel_heights = c(1,0,1), labels = c('A','','C'), ncol=1)
dm_pop_comboplot2
ggsave('Figures/dm-popmodel-plots.pdf', dm_pop_comboplot2, width=4.4, height = 4.2, dpi=300)
ggsave('Figures/dm-popmodel-plots.tiff', dm_pop_comboplot2, width=4.4, height = 4.2, dpi=300)
