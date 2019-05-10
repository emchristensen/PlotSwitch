# Running multiple MARK models for model selection
# Ellen K. Bledsoe
# June 2018
# Modified EMC for Plot Switch project 5/8/19

# LIBRARIES and preliminary parameters
library(dplyr)
library(RMark)

source('FinalAnalysis/population_model_functions.R')
#source('population_model_functions.R')

# first period of capture history desired. the switch happened between periods 436 and 437. the paper begins at period 415
#first_per = 413
first_per = 437

# ==========================================
# read capture history data or create new
# ==========================================

# read in capture history
#mark_trmt_all <- read.csv('Data/capture_history_DM.csv', stringsAsFactors = F)

# create capture history from original data
rdat <- read.csv('PortalData/Rodents/Portal_rodent.csv', header = TRUE, na.strings = c(""), stringsAsFactors = FALSE)
plotswitchplots = c(4,11,14,17,6,13,18,5,7,24)

rdat_filtered = dplyr::filter(rdat, period>=first_per, plot %in% plotswitchplots)
pdat <- data.frame(plot = 1:24, treatment = c('CC','CE','EE','CC','XC','EC',
                                              'XC','CE','CX','XX','CC','CX',
                                              'EC','CC','EE','XX','CC','EC',
                                              'EE','EE','EE','CE','XX','XC'))
#sp = 'PE'
for (sp in c('PM','RO','RM','DM','PP')) {
  
  # filter main data frame by species and run individual tag cleanup
  spdat = individual_tag_cleanup(sp, rdat_filtered)
  spdat_trt = merge(spdat, pdat, by=c('plot'))
  
  # create treatment history file
  mark_trmt_all = create_trmt_hist(spdat_trt)
  
  # ==========================================
  # prep data for RMark
  # ==========================================
  all_ms <- select(mark_trmt_all, ch = captures)
  
  # Process data
  ms.pr = process.data(all_ms, begin.time = first_per, model = "Multistrata")
  
  # Create default design data
  ms.ddl = make.design.data(ms.pr)
  
  # add design covariates for after switch
  after_switch = as.factor(437:476)
  
  ms.ddl$S$after_switch = 0
  ms.ddl$S$after_switch[ms.ddl$S$time %in% after_switch] = 1
  
  ms.ddl$p$after_switch = 0
  ms.ddl$p$after_switch[ms.ddl$p$time %in% after_switch] = 1
  
  ms.ddl$Psi$after_switch = 0
  ms.ddl$Psi$after_switch[ms.ddl$Psi$time %in% after_switch] = 1
  
  # =======================================================
  # Run the models and examine the output
  # =======================================================
  
  #MarkViewer="Notepad" # edit to make results pop up on a Mac
  
  run.ms = function(){
    
    S.dot = list(formula = ~ 1) 
    S.stratum = list(formula =  ~ -1 + stratum)
    S.time = list(formula = ~ -1 + after_switch)
    S.strat_time = list(formula = ~ -1 + stratum + after_switch)
    S.strat_x_time = list(formula = ~ -1 + stratum*after_switch)
    
    p.dot = list(formula = ~ 1) 
    #p.stratum = list(formula =  ~ -1 + stratum)
    #p.time = list(formula = ~ -1 + after_switch)
    #p.strat_time = list(formula = ~ -1 + stratum + after_switch)
    #p.strat_x_time = list(formula = ~ -1 + stratum*after_switch)
    
    Psi.dot = list(formula = ~ 1, link = "logit") 
    Psi.s = list(formula =  ~ -1 + stratum:tostratum, link = "logit")
    Psi.time = list(formula = ~ -1 + after_switch, link = "logit")
    Psi.strat_time = list(formula = ~ -1 + stratum:tostratum + after_switch, link = "logit")
    Psi.strat_x_time = list(formula = ~ -1 + stratum:tostratum*after_switch, link = "logit")
    
    ms.model.list = create.model.list("Multistrata")
    
    ms.results = mark.wrapper(ms.model.list,
                              data = ms.pr, ddl = ms.ddl,
                              options = "SIMANNEAL")
    
    return(ms.results)
    
  }
  
  ms.results = run.ms()
  model.table <- as.data.frame(ms.results$model.table)
  write.csv(model.table, paste0("Data/PopModelSelection_437_476/modelselection_",sp,"_20190510.csv"))
}
#top_model_summary <- ms.results$S.strat_x_time.p.dot.Psi.strat_x_time$results$real
#write.csv(top_model_summary, "top_model_summary_20190508.csv")

