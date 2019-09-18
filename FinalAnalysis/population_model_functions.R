library(igraph)

#' @description used by individual_tag_cleanup: will split apart multiple uses of same tag number if captures are
#'              too far apart in time to be considered the same individual
#' @author Erica Christensen
#' @param indiv
#' @param datframe
#' 
split_dup_tags = function(indiv,datframe,interval) {
  # Function for splitting instances of same tag being used in different time periods
  
  indiv = indiv[order(indiv$date),]           # put in chronological order if not already
  tagrep = 0
  for (n in 2:length(indiv$date)) {
    if (indiv$date[n]-indiv$date[n-1] > interval) {   # if interval between captures > 1.5 years,
      tagrep = tagrep+1
      ids = indiv$Record_ID[n:length(indiv$date)] # find record IDs of subsequent entries to connect back to datframe
      datframe[datframe$Record_ID %in% ids,'tagrep'] = tagrep 
    }
  }                   
  
  assign('datframe',datframe,envir=.GlobalEnv)  # make sure changes to datframe are applied outside this function
}

#' @param sp (character) species of interest
#' @param dat data frame of all rodent data from PortalData
#' @param interval maximum allowable interval between captures before the same tag number should be considered a new individual
#' @param nulls null values
#' 
#' @author Erica Christensen
#' 
individual_tag_cleanup = function(sp, dat, interval=5*365, nulls = c(0,-1,'','000000',9999,NA)) {
  sp_dat = dplyr::filter(dat, species %in% sp)
  
  # remove untagged individuals
  sp_dat = sp_dat[!is.na(sp_dat$tag),]
  
  # add column for full date
  sp_dat$date = as.Date(paste(sp_dat$year,sp_dat$month,sp_dat$day,sep='-'), format='%Y-%m-%d')
  
  #----------------------------------------------------------------------------------------
  # Strip tag and date data from larger data frame
  
  # create data frames with just record_ID, date, and associated tag numbers (tag, ltag, prevrt, prevlet)
  rtag = sp_dat[,c('recordID','date','tag')]
  ltag = sp_dat[,c('recordID','date','ltag')]
  prtag = sp_dat[,c('recordID','date','prevrt')]
  pltag = sp_dat[,c('recordID','date','prevlet')]
  
  # split up previous tags (may contain multiple tags in tag column)
  recordID = c()
  date = c()
  tag = c()
  for (n in 1:length(prtag$prevrt)) {
    rts = unlist(strsplit(prtag$prevrt[n],','))
    for (m in 1:length(rts)) {
      recordID = append(recordID,prtag$recordID[n])
      date = append(date,prtag$date[n])
      tag = append(tag,rts[m])
    }
  }
  prt = data.frame(recordID,date,tag,stringsAsFactors=F)
  
  recordID = c()
  date = c()
  tag = c()
  for (n in 1:length(pltag$prevlet)) {
    lts = unlist(strsplit(pltag$prevlet[n],','))
    for (m in 1:length(lts)) {
      recordID = append(recordID,pltag$recordID[n])
      date = append(date,pltag$date[n])
      tag = append(tag,lts[m])
    }
  }
  plt = data.frame(recordID,date,tag,stringsAsFactors=F)
  
  # remove nulls from these data frames
  rtag = rtag[!rtag$tag %in% nulls,] 
  ltag = ltag[!ltag$ltag %in% nulls,]
  prt = prt[!prt$tag %in% nulls,]
  plt = plt[!plt$tag %in% nulls,]
  
  # Put data frames together into one big data frame
  ltag = dplyr::rename(ltag, tag=ltag)
  datframe = rbind(rtag, ltag, prt, plt)
  datframe$tagrep = 0
  
  # strip leading or trailing whitespace from tags
  datframe$tag <- stringr::str_trim(datframe$tag,side='both')
  
  # make sure all tags are 6 digits (differences in data entry, 001262 versus 1262)
  datframe$tag <- stringr::str_pad(datframe$tag, 6, pad = "0")
  
  #----------------------------------------------------------------------------------------------
  # Split up tags used multiple times in data set
  
  # Loop through unique tag numbers and time-splitting function
  for (t in unique(datframe$tag)) {
    indiv = dplyr::filter(datframe, tag==t) %>% arrange(date)
    if (any(diff(indiv$date)>interval))  {        # if time between captures > 1.5 years, run split_dup_tags
      split_dup_tags(indiv,datframe,interval)
      print(paste('duplicate tag? ',t))
    }
  }
  
  datframe$tagunique = paste(datframe$tag,'_',datframe$tagrep,sep='')
  
  #---------------------------------------------------------------------------------------
  # merge unique tags with original records in data set
  
  R_ID = c()
  utags = c()
  for (id in unique(datframe$recordID)) {
    record = datframe[datframe$recordID==id,]
    R_ID = append(R_ID,id)
    utags = append(utags,paste(record$tagunique,collapse=','))
  }
  
  # Data frame of record IDs and unique tags
  uniquetags = data.frame(recordID=R_ID, tagunique = utags,stringsAsFactors=F)
  
  # Merge with main data frame
  newdatframe = merge(sp_dat,uniquetags,by='recordID')
  
  
  #--------------------------------------------------------------------------------------
  # Isolate unique tag numbers (vertices)
  alltags = unique(datframe$tagunique)
  
  #--------------------------------------------------------------------------------------
  # Isolate connections between tags (L/R/previous) for edges
  
  s = c()   # empty list for storing sets of associated tags
  
  for (n in 1:length(sp_dat$tag)) {        # Loop through all records in data frame
    edge = as.character(unlist(strsplit(newdatframe$tagunique[n],',')))
    edge = edge[!is.na(edge)]          # remove NAs
    edge = edge[!edge %in% nulls]      # remove nulls
    s[n] = list(edge)
  }
  
  #--------------------------------------------------------------------------------------
  # Make it a graph (igraph package)
  
  sp_graph = igraph::graph.empty() + vertices(alltags)
  
  for (edge in s) {
    sp_graph = sp_graph + path(as.vector(edge))
  }
  
  #decomposes into connected components (independent subgraphs) each representing an individual
  graphs = igraph::decompose.graph(sp_graph) 
  
  
  #-----------------------------------------------------------------------------------
  # Group together records associated with each individual (subgraph)
  
  # Make data frame of all tag numbers and which group number (individual) they belong to
  tag = c()
  group = c()
  for (n in 1:length(graphs)) {
    vertices = as.vector(unlist(vertex.attributes(graphs[[n]])))    # vector of tags that represent an individual
    tag = append(tag,vertices)
    group = append(group,rep(n,length(vertices)))
  }
  
  taggroups = data.frame(tagunique=tag,group,stringsAsFactors=F)
  
  # Merge tag/group data frame with capture data
  for (n in 1:length(newdatframe$tagunique)) {
    tgs = as.character(unlist(strsplit(newdatframe$tagunique[n],',')))
    newdatframe$tagunique1[n] = tgs[1]
  }
  new_sp_dat = merge(newdatframe,taggroups,by.x='tagunique1',by.y='tagunique')
  
  # sort by time (record id) and select columns
  selectdat = new_sp_dat %>% 
    dplyr::select(recordID, month, day, year, period, plot, note1, stake, species, sex,
                  age, reprod, testes, vagina, pregnant, nipples, lactation,
                  hfl, wgt, tag, note2, ltag, note3, prevrt, prevlet, note5,
                  date, tagunique, group) %>%
    arrange(recordID)

  return(selectdat)
}


#' @description creates capture history of individual rodents
#' @param dat data frame of rodent captures; after cleaning by individual_tag_cleanup
#' @author Ellen Bledsoe
create_trmt_hist = function(dat) {
  individuals <- unique(dat$group)
  prd <- seq(min(dat$period), max(dat$period))
  
  MARK_data = data.frame("captures" = 1,
                         "censored" = 1,
                         "tags" = 1)
  
  outcount = 0
  
  for (t in 1:length(individuals)) {
    
    capture_history = "" #create empty string
    
    for (p in 1:length(prd)) {
      
      tmp <- dat[which(dat$group == individuals[t] & dat$period == prd[p]),]
      
      if (nrow(tmp) == 0) {
        state = "0"
      } else {
        if (tmp$treatment == 'CC') {
          state = "A"
        } else if (tmp$treatment == 'EC') {
          state = "B"
        } else if (tmp$treatment == 'XC') {
          state = "C"
        } else if (tmp$treatment == 'XX') {
          state = "D"
        } else if (tmp$treatment == 'EE') {
          state = 'E'
        } else if (tmp$treatment == 'CE') {
          state = 'G'
        } else if (tmp$treatment == 'CX') {
          state = 'H'
        }
      }
      capture_history = paste0(capture_history, state)
    }
    
    tmp2 <- which(dat$group == individuals[t])
    censored = 1
    
    outcount = outcount + 1
    MARK_data[outcount, ] <- c(capture_history, censored, individuals[t])
    
  }
  
  return(MARK_data)
  
}


#' @author Ellen Bledsoe
run.ms <- function(ms.pr,
                   S_dot = list(formula = ~ 1),
                   S_stratum = list(formula =  ~ -1 + stratum),
                   p_dot = list(formula =  ~ 1),
                   p_stratum = list(formula =  ~ -1 + stratum),
                   Psi_s = list(formula =  ~ -1 + stratum:tostratum, link = "logit")) {

  # Create default design data
  ms.ddl = make.design.data(ms.pr)

 
  # RMark function for Portal data
  if (is.null(S_dot)) {
    S.stratum = S_stratum
  } else if (is.null(S_stratum)) {
    S.dot = S_dot
  } else {
    S.stratum = S_stratum
    S.dot = S_dot
  }

  if (is.null(p_dot)) {
    p.stratum = p_stratum
  } else if (is.null(p_stratum)) {
    p.dot = p_dot
  } else {
    p.stratum = p_stratum
    p.dot = p_dot
  }

  Psi.s = Psi_s

  # Create model list and run assortment of models
  ms.model.list = create.model.list("Multistrata")

  ms.results = mark.wrapper(ms.model.list,
                            data = ms.pr, ddl = ms.ddl,
                            options="SIMANNEAL")

  # Return model table and list of models
  return(ms.results)

}



#' @title run pop model
#' @author Erica Christensen
#' @description this is a wrapper function that runs tag cleanup, the chosen rmark model, and 
#'              writes the output to csv for later plotting/analysis
#' @param rdat data frame of rodent data, filtered to the appropriate time period and plots
#' @param sp species (2-letter character code)
#' @param write_cap_history T/F whether to write capture history to csv
#' @param date_run string: date of the run for ID purposes. will be put in input/output file names.
run_species_pop_model <- function(rdat, sp, write_cap_history = F, date_run) {
  # plot treatment information (create data frame)
  pdat <- data.frame(plot = 1:24, treatment = c('CC','CE','EE','CC','XC','EC',
                                                'XC','CE','CX','XX','CC','CX',
                                                'EC','CC','EE','XX','CC','EC',
                                                'EE','EE','EE','CE','XX','XC'))
  # filter main data frame by species and run individual tag cleanup
  spdat = individual_tag_cleanup(sp, rdat)
  spdat_trt = merge(spdat, pdat, by=c('plot'))
  
  # number of individuals
  n_indiv = length(unique(spdat$group))
  
  # create treatment history file
  mark_sp = create_trmt_hist(spdat_trt)
  # resulting warnings are probably for animals captured twice in same sampling period
  
  # write to csv
  if (write_cap_history == T) {
    write.csv(mark_sp, file=paste0('Data/capture_history_',sp,'.csv'), row.names = F)
  }
  
  # read info on best model from model selection file
  bestmodel <- read.csv(paste0('Data/PopModelSelection_afteronly/modelselection_',sp,'_',date_run,'.csv'),stringsAsFactors = F, nrows=1)
  
  # prep data for RMark
  all_ms <- select(mark_sp, ch = captures)
  first_per <- min(spdat_trt$period)
  
  # Process data
  ms.pr = process.data(all_ms, begin.time = first_per, model = "Multistrata")
  
  # run model just on after-switch data
  ms.results = run.ms(ms.pr,
                      S_dot = NULL,
                      S_stratum = list(formula = as.formula(bestmodel$S)),
                      p_dot = list(formula = ~ 1),
                      p_stratum = NULL,
                      Psi_s = list(formula = as.formula(bestmodel$Psi), link = "logit"))
  
  # write output to csv
  rmark_results <- ms.results$S.stratum.p.dot.Psi.s$results$real
  rmark_results$n_indiv <- n_indiv
  write.csv(rmark_results, paste0("Data/PopModelBest_afteronly/MARK_", sp, "_top_model_summary_",date_run,".csv"))
}

#' @author Ellen Bledsoe
prep_RMark_data_for_plotting <- function(data){

  # add descriptive columns
  data$metric = gsub(' .*$', '', data$X)
  data$stratum = gsub(' g1.*$', '', data$X)
  data$stratum = gsub('.* s', '', data$stratum)

  data$Treatment = gsub('C', ' Rodent+', data$stratum)
  data$Treatment = gsub('A', ' Control', data$Treatment)
  data$Treatment = gsub('B', ' Kangaroo rat+', data$Treatment)
  data$Treatment = gsub('to', '\nto\n', data$Treatment)
  
  outdata <- data %>% 
    filter(metric %in% c("S","Psi"))
  
  return(outdata)
  
}

#' @author Ellen Bledsoe, modified by Erica Christensen
#' @param data
#' @param colorvalues color values for manual color scale (make sure you have the right number)
#' @param maintitle title for plot
plot_estimated_survival <- function(data, colorvalues, maintitle){
  
  # plot estimated survival metrics from RMark

  if (data$metric[1] == 'S') {
    y_label = "Estimated survival"
    x_var = 'Treatment'
  } else if (data$metric[1] == 'Psi') {
    y_label = "Transition probability"
    x_var = 'Transition'
  } else {
    y_label= ''
  }
  
  plot <- ggplot(data, color = Treatment) +
    geom_pointrange(aes(x = eval(parse(text=x_var)), y = estimate, 
                        ymin = (estimate - se), ymax = (estimate + se), 
                        color = Treatment), 
                    position = position_dodge(.1), size = .5) +
    scale_colour_manual(name = 'Treatment:',
                        values = colorvalues) + 
    ggtitle(maintitle) +
    ylab(y_label) + 
    xlab('') +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 6),
          axis.title.y = element_text(size = 10),
          #axis.text.x = element_text(angle=30),
          legend.position = "none",
          plot.margin = margin(l = 5, t = 0, b = 0, r = 5))
  
  return(plot)
  
}



#' @title new captures of a given species by year since treatment change
#' @author Erica Christensen
#' 
#' @param species (2 letter species code)
#' @param rdat original unfiltered rdat dataframe
#' @param tdat trapping data from Portal_rodent_trapping.csv
new_captures_by_year <- function(sp, rdat, tdat) {
  # create df of individuals going back to 2012 to make sure we see the first capture of all animals seen in 2015-2018. Includes all plots
  sp_only = individual_tag_cleanup(sp, dplyr::filter(rdat, period>= 403))
  
  # create empty data frame
  first_captures <- setNames(data.frame(matrix(ncol=length(names(sp_only)), nrow=0)), names(sp_only))
  
  # loop through individuals and identify first capture
  for (g in unique(sp_only$group)) {
    indiv <- dplyr::filter(sp_only, group==g)
    first <- dplyr::filter(indiv, period==min(indiv$period)) %>% dplyr::slice(1)
    first_captures <- rbind(first_captures, first)
  }
  
  # get tally of new individuals by period and plot, and filter for post-switch periods and relevant plots
  new_per_plot <- first_captures %>%
    dplyr::select(period, date, plot) %>%
    group_by(period, date, plot) %>%
    summarise(count = n()) %>%
    arrange(period, plot) %>%
    dplyr::filter(period>=437, plot %in% plotswitchplots) %>%
    ungroup()
  
  alldateplots = new_per_plot %>% expand(period,plot) %>% merge(new_per_plot, all.x=T) %>% merge(tdat, by=c('period','plot'), all.x=T)
  alldateplots$count[(is.na(alldateplots$count) & alldateplots$sampled==1)] <- 0
  
  # merge with plot treatment df and get avg new individuals
  new_per_plot_trt = merge(alldateplots, pdat, all.x=T) %>% dplyr::select(plot, period, count, day, month, year, sampled, treatment)
  
  # remove rows where plots were not sampled
  new_per_plot_trt = new_per_plot_trt[new_per_plot_trt$sampled==1,]
  
  # calculate average number of new per year after switch
  new_per_plot_trt$yr = new_per_plot_trt$year
  new_per_plot_trt$yr[new_per_plot_trt$month<4] <- new_per_plot_trt$yr[new_per_plot_trt$month<4]-1
  
  # create column for "years after switch"
  new_per_plot_trt$yr_since_change <- NA
  new_per_plot_trt$yr_since_change[new_per_plot_trt$period %in% 437:447] <- 1
  new_per_plot_trt$yr_since_change[new_per_plot_trt$period %in% 448:460] <- 2
  new_per_plot_trt$yr_since_change[new_per_plot_trt$period %in% 461:472] <- 3
  new_per_plot_trt$yr_since_change[new_per_plot_trt$period > 472] <- 4 # April-July 2018
  
  # create column to arrange treatments on xaxis for plotting
  xaxisarrangement = data.frame(yr_since_change = c(1,1,1,2,2,2,3,3,3,4,4,4),
                               treatment = rep(c('CC','EC','XC'),4),
                               xposition = c(1,1.2,1.4,2,2.2,2.4,3,3.2,3.4,4,4.2,4.4))
  #new_per_plot_trt2 = merge(new_per_plot_trt, xaxisarrangement)
  
  # calculate total new animals per treatment type in each year since switch
  t <- new_per_plot_trt %>% group_by(yr_since_change, treatment, plot) %>%
    summarise(total_new_per_year=sum(count)) %>%
    merge(xaxisarrangement)  %>%
    dplyr::filter(yr_since_change<4) %>%
    group_by(yr_since_change, treatment) %>%
    mutate(mean_by_trt = mean(total_new_per_year))
  
  # NOTE: there are 4 CC plots and only 3 of each EC and XC plots
  
  return(t)
}
