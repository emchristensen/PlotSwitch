# Script for cleaning up individual tag data from Portal data set
# Erica Christensen 11/2014
# update 2/2016
# update 4/2019 for use with PlotSwitch project



rm(list=ls(all=TRUE))

#library(stringr)
library(igraph)
library(dplyr)

###########################################################################################
# read in data
###########################################################################################

dat = read.csv('PortalData/Rodents/Portal_rodent.csv',as.is=TRUE,colClasses=c(note1='character',
                                                                         tag='character'),na.strings = '')


# remove negative period codes -- not part of regular census
dat = dplyr::filter(dat, period>0)
# restrict data to time period and plots of interest
selectedplots = c(4, 5, 6, 7, 11, 13, 14, 17, 18, 24)
dat = dplyr::filter(dat, period %in% 415:476,
                    plot %in% selectedplots)

#############################################################################################
# some parameters and functions for use later
#############################################################################################
interval = 5*365      # maximum acceptable time (in days) for an animal to disappear between captures, 
                        # for spearating multiple uses of same tag number
lifespan = 5*365        # maximum life span for selected species
# declare nulls
nulls = c(0,-1,'','000000',9999,NA)

split_dup_tags = function(indiv,datframe) {
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

############################################################################################
# main function
############################################################################################
individual_tag_cleanup = function(sp) {
  sp_dat = dplyr::filter(dat, species==sp)
  
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
  
  # remove nulls from these data frames
  rtag = rtag[!rtag$tag %in% nulls,] 
  ltag = ltag[!ltag$ltag %in% nulls,]
  prtag = prtag[!prtag$prevrt %in% nulls,]
  pltag = pltag[!pltag$prevlet %in% nulls,]
  
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
      split_dup_tags(indiv,datframe)
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

  sp_graph = graph.empty() + vertices(alltags)
  
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
  
  # sort by time (record id)
  new_sp_dat = new_sp_dat[order(new_sp_dat$recordID),]
  

  
  selectdat = new_sp_dat[,c('recordID','month','day','year','period','plot','species','sex','age',
                            'reprod','testes','vagina','pregnant','nipples','lactation','hfl',
                            'wgt','tag','note2','ltag','note3','prevrt','prevlet','note5','date',
                            'tagunique','group')]
  
  return(selectdat)
}

###################################################################
# Loop through species
####################################################################
for (species in c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PM','PP','RM','SF','SH')) {
  selectdat = individual_tag_cleanup(species)
  write.csv(selectdat,paste('Data/', sp,'_indiv.csv',sep=''),row.names=F)
}


#---------------------------------------------------------------------------------
# Look for errors in grouping

# multiple captures on same day (probably 2 or more individuals)
tag_check = c()
for (grp in unique(new_sp_dat$group)) {
  indiv = new_sp_dat[new_sp_dat$group==grp,]
  if (length(indiv$date) != length(unique(indiv$date))) {
    tag_check = rbind(tag_check,grp)
  }
}

# Captures after death
for (grp in unique(new_sp_dat$group)) {
  indiv = new_sp_dat[new_sp_dat$group==grp,]
  if ('D' %in% indiv$note5) {
    death = indiv[indiv$note5=='D',]
    if (any(indiv$date > max(death$date,na.rm=T))) {
      tag_check = rbind(tag_check,grp)
    }
  }
}

# Conflicting sexual characteristics
for (grp in unique(new_sp_dat$group)) {
  indiv = new_sp_dat[new_sp_dat$group==grp,]
  m_char = (('R' %in% indiv$testes) || ('S' %in% indiv$testes) || ('M' %in% indiv$testes))
  f_char = (('S' %in% indiv$vagina) || ('P' %in% indiv$vagina) || ('P' %in% indiv$pregnant) || 
              ('R' %in% indiv$nipples) || ('E' %in% indiv$nipples) || ('B' %in% indiv$nipples) || 
              ('L' %in% indiv$lactation))
  if (f_char && m_char) {
    tag_check = rbind(tag_check,grp)
  }
}

# capture history exceeds expected lifespan of animal
for (grp in unique(new_sp_dat$group)) {
  indiv = new_sp_dat[new_sp_dat$group==grp,]
  if (max(indiv$date)-min(indiv$date) > lifespan) {
    tag_check = rbind(tag_check,grp)
  }
}
# ----------------------------------------------------------
# CSV output

# groups to remove from dataset
remv = unique(tag_check)
for (n in remv) {
  new_sp_dat[new_sp_dat$group==n,'group'] = -n
}
