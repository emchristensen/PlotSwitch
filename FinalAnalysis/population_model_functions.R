library(igraph)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#' @description used by individual_tag_cleanup: will split apart multiple uses of same tag number if captures are
#'              too far apart in time to be considered the same individual
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
    dplyr::select(recordID, month, day, year, period, plot, species, sex,
                  age, reprod, testes, vagina, pregnant, nipples, lactation,
                  hfl, wgt, tag, note2, ltag, note3, prevrt, prevlet, note5,
                  date, tagunique, group) %>%
    arrange(recordID)

  return(selectdat)
}


#' @description creates capture history of individual rodents
#' @param dat data frame of rodent captures; after cleaning by individual_tag_cleanup
#' Written by Ellen Bledsoe
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


#' written by Ellen Bledsoe
run.ms <- function(S_dot = list(formula = ~ 1), 
                   S_stratum = list(formula =  ~ -1 + stratum + after_switch), 
                   p_dot = list(formula =  ~ 1), 
                   p_stratum = list(formula =  ~ -1 + stratum + after_switch), 
                   Psi_s = list(formula =  ~ -1 + stratum:tostratum + after_switch, link = "logit")) {
  
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

#' written by Ellen Bledsoe
prep_RMark_data_for_plotting <- function(data){

  # add descriptive columns
  data$metric = gsub(' .*$', '', data$X)
  data$stratum = gsub(' g1.*$', '', data$X)
  data$stratum = gsub('.* s', '', data$stratum)

  data$Treatment = gsub('C', 'Former rodent \nremoval', data$stratum)
  data$Treatment = gsub('A', 'Long-term \ncontrol', data$Treatment)
  data$Treatment = gsub('B', 'Former kangaroo \nrat removal', data$Treatment)
  
  outdata <- data %>% 
    filter(metric == "S")
  
  return(outdata)
  
}

# written by Ellen Bledsoe
plot_estimated_survival <- function(data, maintitle){
  
  # plot estimated survival metrics from RMark
  
  #x_axis_title <- expression(paste(italic("C. baileyi"), " establishment"))
  
  plot <- ggplot(data, color = Treatment) +
    geom_pointrange(aes(x = Treatment, y = estimate, 
                        ymin = (estimate - se), ymax = (estimate + se), 
                        color = Treatment), 
                    position = position_dodge(.1), size = .75) +
    scale_colour_manual(values = cbbPalette) + 
    ggtitle(maintitle) +
    ylab("Estimated survival") + 
    theme_classic() +
    theme(panel.border = element_rect(fill = NA, colour = "black", size = 1.25),
          plot.subtitle = element_text(size = 14, hjust = -.4),
          axis.line = element_line(size = .25),
          axis.title.x = element_text(size = 12, margin = margin(r = 10)),
          axis.title.y = element_text(size = 12, margin = margin(r = 10)),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "none",
          legend.title = element_blank(),
          plot.margin = margin(l = 5, t = 20))
  
  return(plot)
  
}
